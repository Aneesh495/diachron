//! Native Rust MCMC chain over a timed tree and binary gain/loss parameters.
//!
//! The Python layer drives the research workflow, but the inner sampling loop is
//! available natively here so a chain can run without crossing the language
//! boundary on every step. The chain samples the tree topology and node times
//! together with the gain and loss rates against a target that is the character
//! pruning likelihood plus a coalescent tree prior plus exponential rate priors.
//! It uses the same serializable RNG as the Python driver, so a native chain is
//! equally reproducible and checkpointable.
//!
//! This module is intentionally self-contained: it owns its tree copy, applies
//! moves in place with revert-on-reject, and exposes thinned scalar draws. It is
//! validated against the pruning oracle (the likelihood it maximizes is the same
//! one the enumerator checks) and against the Python chain's behavior on shared
//! seeds for the scalar moves.

use crate::pruning::log_likelihood_block;
use crate::rng::Rng;
use crate::subst::BinaryGainLoss;
use crate::tree::Tree;

/// Target density pieces and data needed to score a state.
pub struct Target {
    /// Character columns in taxon order; negative entries are missing.
    pub columns: Vec<Vec<Option<usize>>>,
    /// Coalescent effective scale (years), for the tree prior.
    pub pop_size: f64,
}

/// One native chain.
pub struct NativeChain {
    pub tree: Tree,
    pub gain: f64,
    pub loss: f64,
    pub rng: Rng,
    target: Target,
    log_post: f64,
    scale_lambda: f64,
    time_window_frac: f64,
    accepts: u64,
    proposals: u64,
    pub draws: Vec<ScalarDraw>,
}

#[derive(Clone, Debug)]
pub struct ScalarDraw {
    pub iteration: u64,
    pub log_post: f64,
    pub root_age: f64,
    pub gain: f64,
    pub loss: f64,
}

impl NativeChain {
    pub fn new(tree: Tree, gain: f64, loss: f64, seed: u64, target: Target) -> Self {
        let mut chain = NativeChain {
            tree,
            gain,
            loss,
            rng: Rng::seed_from_u64(seed),
            target,
            log_post: f64::NEG_INFINITY,
            scale_lambda: 0.5,
            time_window_frac: 0.1,
            accepts: 0,
            proposals: 0,
            draws: Vec::new(),
        };
        chain.log_post = chain.evaluate();
        chain
    }

    fn root_prior(&self) -> [f64; 2] {
        let z = self.gain + self.loss;
        [self.loss / z, self.gain / z]
    }

    /// Unnormalized log posterior of the current state.
    fn evaluate(&self) -> f64 {
        if self.gain <= 0.0 || self.loss <= 0.0 {
            return f64::NEG_INFINITY;
        }
        if !self.tree.times_valid() {
            return f64::NEG_INFINITY;
        }
        let model = BinaryGainLoss::new(self.gain, self.loss);
        let prior = self.root_prior();
        let ll = log_likelihood_block(&self.tree, &model, &self.target.columns, &prior);
        if !ll.is_finite() {
            return f64::NEG_INFINITY;
        }
        ll + self.coalescent_log_prior() + self.rate_log_prior()
    }

    fn rate_log_prior(&self) -> f64 {
        // Exponential(mean 1e-3) on each rate.
        let rate: f64 = 1.0 / 1e-3;
        rate.ln() - rate * self.gain + rate.ln() - rate * self.loss
    }

    fn coalescent_log_prior(&self) -> f64 {
        // Constant-size coalescent on internal node times.
        let n = self.tree.nodes.len();
        let mut is_parent = vec![false; n];
        for node in &self.tree.nodes {
            if let Some(p) = node.parent {
                is_parent[p] = true;
            }
        }
        let mut events: Vec<(f64, bool)> = (0..n)
            .map(|i| (self.tree.nodes[i].height, is_parent[i]))
            .collect();
        events.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut lp = 0.0;
        let mut lineages: i64 = 0;
        let mut prev_t = 0.0;
        for (t, is_coal) in events {
            let dt = t - prev_t;
            if lineages >= 2 && dt > 0.0 {
                let rate = (lineages * (lineages - 1)) as f64 / (2.0 * self.target.pop_size);
                lp += -rate * dt;
            }
            if is_coal {
                if lineages >= 2 {
                    let rate =
                        (lineages * (lineages - 1)) as f64 / (2.0 * self.target.pop_size);
                    lp += rate.ln();
                }
                lineages -= 1;
            } else {
                lineages += 1;
            }
            prev_t = t;
        }
        lp
    }

    fn accept(&mut self, log_post_new: f64, log_hastings: f64) -> bool {
        if !log_post_new.is_finite() {
            return false;
        }
        let log_alpha = log_post_new - self.log_post + log_hastings;
        if log_alpha >= 0.0 {
            return true;
        }
        self.rng.next_f64().max(1e-300).ln() < log_alpha
    }

    /// One Metropolis-within-Gibbs step.
    pub fn step(&mut self) {
        let kind = self.rng.next_below(5);
        self.proposals += 1;
        let mut accepted = false;
        match kind {
            0 => {
                // scale gain
                let u = self.rng.next_f64();
                let factor = (self.scale_lambda * (u - 0.5)).exp();
                let old = self.gain;
                self.gain *= factor;
                let lp = self.evaluate();
                if self.gain > 0.0 && self.accept(lp, factor.ln()) {
                    self.log_post = lp;
                    accepted = true;
                } else {
                    self.gain = old;
                }
            }
            1 => {
                // scale loss
                let u = self.rng.next_f64();
                let factor = (self.scale_lambda * (u - 0.5)).exp();
                let old = self.loss;
                self.loss *= factor;
                let lp = self.evaluate();
                if self.loss > 0.0 && self.accept(lp, factor.ln()) {
                    self.log_post = lp;
                    accepted = true;
                } else {
                    self.loss = old;
                }
            }
            2 => accepted = self.slide_internal_node(),
            3 => accepted = self.scale_root(),
            _ => accepted = self.nni(),
        }
        if accepted {
            self.accepts += 1;
        }
    }

    fn internal_nodes(&self) -> Vec<usize> {
        (0..self.tree.nodes.len())
            .filter(|&i| !self.tree.is_tip(i))
            .collect()
    }

    fn slide_internal_node(&mut self) -> bool {
        let internals = self.internal_nodes();
        if internals.is_empty() {
            return false;
        }
        let node = internals[self.rng.next_below(internals.len())];
        let lower = self.tree.nodes[node]
            .children
            .iter()
            .map(|&c| self.tree.nodes[c].height)
            .fold(f64::NEG_INFINITY, f64::max);
        let upper = match self.tree.nodes[node].parent {
            Some(p) => self.tree.nodes[p].height,
            None => f64::INFINITY,
        };
        let window = self.time_window_frac * self.tree.nodes[self.tree.root].height.max(1.0);
        let step = (self.rng.next_f64() - 0.5) * 2.0 * window;
        let proposed = self.tree.nodes[node].height + step;
        if proposed <= lower + 1e-9 || (upper.is_finite() && proposed >= upper - 1e-9) {
            return false;
        }
        let old = self.tree.nodes[node].height;
        self.tree.nodes[node].height = proposed;
        self.tree.set_branch_lengths_from_heights();
        let lp = self.evaluate();
        if self.accept(lp, 0.0) {
            self.log_post = lp;
            true
        } else {
            self.tree.nodes[node].height = old;
            self.tree.set_branch_lengths_from_heights();
            false
        }
    }

    fn scale_root(&mut self) -> bool {
        let internals = self.internal_nodes();
        let u = self.rng.next_f64();
        let factor = (self.scale_lambda * (u - 0.5)).exp();
        let saved: Vec<f64> = internals.iter().map(|&i| self.tree.nodes[i].height).collect();
        for &i in &internals {
            self.tree.nodes[i].height *= factor;
        }
        self.tree.set_branch_lengths_from_heights();
        if !self.tree.times_valid() {
            for (k, &i) in internals.iter().enumerate() {
                self.tree.nodes[i].height = saved[k];
            }
            self.tree.set_branch_lengths_from_heights();
            return false;
        }
        let lp = self.evaluate();
        if self.accept(lp, internals.len() as f64 * factor.ln()) {
            self.log_post = lp;
            true
        } else {
            for (k, &i) in internals.iter().enumerate() {
                self.tree.nodes[i].height = saved[k];
            }
            self.tree.set_branch_lengths_from_heights();
            false
        }
    }

    fn nni(&mut self) -> bool {
        // Collect (parent, internal-child) candidates.
        let n = self.tree.nodes.len();
        let mut candidates = Vec::new();
        for c in 0..n {
            if self.tree.is_tip(c) {
                continue;
            }
            if let Some(p) = self.tree.nodes[c].parent {
                if self.tree.nodes[p].children.len() >= 2 {
                    candidates.push((p, c));
                }
            }
        }
        if candidates.is_empty() {
            return false;
        }
        let (p, c) = candidates[self.rng.next_below(candidates.len())];
        let siblings: Vec<usize> = self.tree.nodes[p]
            .children
            .iter()
            .copied()
            .filter(|&x| x != c)
            .collect();
        if siblings.is_empty() {
            return false;
        }
        let sib = siblings[self.rng.next_below(siblings.len())];
        if self.tree.nodes[c].children.is_empty() {
            return false;
        }
        let gc = self.tree.nodes[c].children[self.rng.next_below(self.tree.nodes[c].children.len())];
        if self.tree.nodes[c].height <= self.tree.nodes[sib].height + 1e-9 {
            return false;
        }
        if self.tree.nodes[p].height <= self.tree.nodes[gc].height + 1e-9 {
            return false;
        }

        // Apply swap.
        self.apply_nni(p, c, sib, gc);
        let lp = self.evaluate();
        if self.accept(lp, 0.0) {
            self.log_post = lp;
            true
        } else {
            // Inverse swap restores the tree.
            self.apply_nni(p, c, gc, sib);
            false
        }
    }

    fn apply_nni(&mut self, p: usize, c: usize, sib: usize, gc: usize) {
        if let Some(pos) = self.tree.nodes[p].children.iter().position(|&x| x == sib) {
            self.tree.nodes[p].children[pos] = gc;
        }
        if let Some(pos) = self.tree.nodes[c].children.iter().position(|&x| x == gc) {
            self.tree.nodes[c].children[pos] = sib;
        }
        self.tree.nodes[sib].parent = Some(c);
        self.tree.nodes[gc].parent = Some(p);
        self.tree.set_branch_lengths_from_heights();
    }

    /// Run warmup (with adaptation) then sampling, recording thinned draws.
    pub fn run(&mut self, warmup: u64, iterations: u64, thin: u64) {
        for it in 0..(warmup + iterations) {
            self.step();
            if it < warmup && self.proposals % 100 == 0 {
                self.adapt();
            }
            if it >= warmup && (it - warmup) % thin == 0 {
                self.draws.push(ScalarDraw {
                    iteration: it,
                    log_post: self.log_post,
                    root_age: self.tree.nodes[self.tree.root].height,
                    gain: self.gain,
                    loss: self.loss,
                });
            }
        }
    }

    fn adapt(&mut self) {
        let rate = self.accepts as f64 / self.proposals.max(1) as f64;
        let lr = 0.5 / (1.0 + self.proposals as f64 / 1000.0);
        let adjust = (lr * (rate - 0.3)).exp();
        self.scale_lambda = (self.scale_lambda * adjust).clamp(1e-3, 5.0);
        self.time_window_frac = (self.time_window_frac * adjust).clamp(1e-4, 0.9);
    }

    pub fn acceptance_rate(&self) -> f64 {
        self.accepts as f64 / self.proposals.max(1) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::Tree;

    fn star_tree(n_tips: usize, root_age: f64) -> Tree {
        // A simple caterpillar so there is topology to explore.
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        tree.nodes[root].height = root_age;
        let mut current = root;
        let mut age = root_age;
        for t in 0..n_tips.saturating_sub(1) {
            age *= 0.7;
            let internal = tree.add_node(Some(current), 0.0, None);
            tree.nodes[internal].height = age;
            tree.add_node(Some(current), 0.0, Some(t));
            current = internal;
        }
        tree.nodes[current].taxon = Some(n_tips - 1);
        tree.root = root;
        tree.set_branch_lengths_from_heights();
        tree
    }

    #[test]
    fn native_chain_runs_and_mixes() {
        let n = 6;
        let tree = star_tree(n, 5000.0);
        // Two informative characters.
        let columns = vec![
            (0..n).map(|i| Some((i % 2) as usize)).collect::<Vec<_>>(),
            (0..n).map(|i| Some(((i / 2) % 2) as usize)).collect::<Vec<_>>(),
        ];
        let target = Target { columns, pop_size: 3000.0 };
        let mut chain = NativeChain::new(tree, 5e-4, 1e-3, 42, target);
        chain.run(500, 1500, 10);
        assert!(!chain.draws.is_empty());
        assert!(chain.acceptance_rate() > 0.0 && chain.acceptance_rate() < 1.0);
        // Log posterior should be finite throughout the retained draws.
        assert!(chain.draws.iter().all(|d| d.log_post.is_finite()));
    }

    #[test]
    fn native_chain_reproducible() {
        let build = || {
            let n = 5;
            let tree = star_tree(n, 4000.0);
            let columns = vec![(0..n).map(|i| Some((i % 2) as usize)).collect::<Vec<_>>()];
            let target = Target { columns, pop_size: 2000.0 };
            NativeChain::new(tree, 5e-4, 1e-3, 7, target)
        };
        let mut a = build();
        let mut b = build();
        a.run(200, 400, 5);
        b.run(200, 400, 5);
        assert_eq!(a.draws.len(), b.draws.len());
        for (da, db) in a.draws.iter().zip(b.draws.iter()) {
            assert!((da.log_post - db.log_post).abs() < 1e-12);
            assert!((da.gain - db.gain).abs() < 1e-15);
        }
    }
}
