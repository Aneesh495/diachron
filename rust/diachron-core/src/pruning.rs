//! Felsenstein pruning likelihood with log-space scaling.
//!
//! Given a rooted tree, a substitution model, and observed states at the tips,
//! `log_likelihood` computes the probability of the observed character under the
//! model, marginalizing over all internal node states. The recursion stores, for
//! each node and each state s, the conditional likelihood
//!
//!     L_node(s) = P(data below node | node is in state s),
//!
//! computed bottom up. For a tip the vector is an indicator of the observed
//! state; a missing observation is the all-ones vector (every state is equally
//! consistent with "we did not observe this cell"), which is what makes missing
//! and absent produce different likelihoods. At an internal node with children
//! c, the per-child contribution is sum_j P_ij(t_c) L_c(j), and the node vector
//! is the product of child contributions.
//!
//! Numerical stability: products of many small probabilities underflow. At each
//! internal node we factor out the maximum entry of the conditional vector into
//! an accumulated log-scale, so the stored vector stays O(1). The final log
//! likelihood adds the accumulated scales back.
//!
//! Multifurcations are handled directly: an internal node simply has more than
//! two children and the product runs over all of them.

use crate::linalg::Mat;
use crate::subst::SubstModel;
use crate::tree::Tree;

/// Observed tip states. `states[taxon]` is `Some(s)` for an observed state `s`
/// in `0..n_states`, or `None` for a missing observation.
pub struct TipStates<'a> {
    pub states: &'a [Option<usize>],
}

/// Compute the log-likelihood of a single character on a tree under a model.
///
/// `root_prior` is the probability distribution over the root state (usually the
/// model's stationary distribution). The returned value is a natural log.
pub fn log_likelihood(
    tree: &Tree,
    model: &dyn SubstModel,
    tips: &TipStates,
    root_prior: &[f64],
) -> f64 {
    PruningPlan::new(tree, model, 1.0).score(tips.states, root_prior)
}

/// Tree traversal and transitions shared across sites of one rate category.
pub(crate) struct PruningPlan<'a> {
    pub(crate) tree: &'a Tree,
    pub(crate) order: Vec<usize>,
    pub(crate) transitions: Vec<Mat>,
    n_states: usize,
}

impl<'a> PruningPlan<'a> {
    pub(crate) fn new(tree: &'a Tree, model: &dyn SubstModel, rate: f64) -> Self {
        Self {
            tree,
            order: tree.postorder(),
            transitions: tree
                .nodes
                .iter()
                .map(|node| model.transition(node.branch_len * rate))
                .collect(),
            n_states: model.n_states(),
        }
    }

    pub(crate) fn score(&self, states: &[Option<usize>], root_prior: &[f64]) -> f64 {
        let tree = self.tree;
        let n_states = self.n_states;
        assert_eq!(root_prior.len(), n_states);
        let n_nodes = tree.nodes.len();
        // Conditional likelihood vectors, one length-n_states block per node.
        let mut cond = vec![0.0f64; n_nodes * n_states];
        // Accumulated log scale per node.
        let mut logscale = vec![0.0f64; n_nodes];

        for &node in &self.order {
            let base = node * n_states;
            if tree.is_tip(node) {
                let taxon = tree.nodes[node].taxon.expect("tip without taxon");
                match states[taxon] {
                    Some(s) => {
                        for k in 0..n_states {
                            cond[base + k] = if k == s { 1.0 } else { 0.0 };
                        }
                    }
                    None => {
                        // Missing: consistent with every state.
                        for k in 0..n_states {
                            cond[base + k] = 1.0;
                        }
                    }
                }
                logscale[node] = 0.0;
            } else {
                // Start from all-ones and multiply in each child's contribution.
                for k in 0..n_states {
                    cond[base + k] = 1.0;
                }
                let mut child_scale = 0.0;
                for &child in &tree.nodes[node].children {
                    let p = &self.transitions[child];
                    let cbase = child * n_states;
                    // partial[i] = sum_j P_ij(t) * L_child(j)
                    for i in 0..n_states {
                        let mut acc = 0.0;
                        for j in 0..n_states {
                            acc += p.at(i, j) * cond[cbase + j];
                        }
                        cond[base + i] *= acc;
                    }
                    child_scale += logscale[child];
                    // Scale after each child so a wide polytomy stays finite too.
                    let block = &mut cond[base..base + n_states];
                    let maxv = block.iter().copied().fold(0.0, f64::max);
                    if maxv <= 0.0 {
                        return f64::NEG_INFINITY;
                    }
                    for value in block {
                        *value /= maxv;
                    }
                    child_scale += maxv.ln();
                }
                logscale[node] = child_scale;
            }
        }

        // Combine at the root with the prior.
        let rbase = tree.root * n_states;
        let mut site = 0.0;
        for k in 0..n_states {
            site += root_prior[k] * cond[rbase + k];
        }
        if site <= 0.0 {
            return f64::NEG_INFINITY;
        }
        site.ln() + logscale[tree.root]
    }
}

/// Log-likelihood of a block of characters that share the same tree and model,
/// summed over sites. `columns[c][taxon]` is the state of taxon at character c.
pub fn log_likelihood_block(
    tree: &Tree,
    model: &dyn SubstModel,
    columns: &[Vec<Option<usize>>],
    root_prior: &[f64],
) -> f64 {
    let plan = PruningPlan::new(tree, model, 1.0);
    columns.iter().map(|col| plan.score(col, root_prior)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subst::BinaryGainLoss;

    /// Build a two-tip tree: root -> {tipA, tipB}, each branch length t.
    fn two_tip(t: f64) -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        tree.add_node(Some(root), t, Some(0));
        tree.add_node(Some(root), t, Some(1));
        tree.root = root;
        tree
    }

    #[test]
    fn two_tip_matches_hand_computation() {
        // With both tips present (state 1), the likelihood is
        // sum_r pi_r P_r1(t) P_r1(t).
        let model = BinaryGainLoss::new(0.4, 0.8);
        let t = 0.6;
        let tree = two_tip(t);
        let pi = model.stationary();
        let p = model.transition(t);
        let expected = pi[0] * p.at(0, 1) * p.at(0, 1) + pi[1] * p.at(1, 1) * p.at(1, 1);
        let states = vec![Some(1usize), Some(1usize)];
        let tips = TipStates { states: &states };
        let got = log_likelihood(&tree, &model, &tips, &pi);
        assert!(
            (got - expected.ln()).abs() < 1e-10,
            "{got} vs {}",
            expected.ln()
        );
    }

    #[test]
    fn missing_differs_from_absent() {
        let model = BinaryGainLoss::new(0.4, 0.8);
        let t = 0.6;
        let tree = two_tip(t);
        let pi = model.stationary();

        let absent = vec![Some(1usize), Some(0usize)];
        let missing = vec![Some(1usize), None];
        let ll_absent = log_likelihood(&tree, &model, &TipStates { states: &absent }, &pi);
        let ll_missing = log_likelihood(&tree, &model, &TipStates { states: &missing }, &pi);
        // Missing marginalizes over tip B, so its likelihood is strictly larger
        // than pinning B to a specific state.
        assert!(ll_missing > ll_absent, "{ll_missing} !> {ll_absent}");
    }

    #[test]
    fn scaling_stable_on_wide_polytomy() {
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let n = 2000;
        for taxon in 0..n {
            tree.add_node(Some(root), 100.0, Some(taxon));
        }
        let model = BinaryGainLoss::new(0.5, 0.5);
        let columns = vec![vec![Some(1); n], vec![Some(0); n]];
        let expected = -((2 * n) as f64) * 2.0f64.ln();
        let score = log_likelihood_block(&tree, &model, &columns, &[0.5, 0.5]);
        assert!((score - expected).abs() < 1e-8, "{score} vs {expected}");
    }

    #[test]
    fn scaling_stable_on_deep_tree() {
        // A caterpillar tree deep enough to underflow without scaling.
        let model = BinaryGainLoss::new(0.5, 0.5);
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let mut current = root;
        let mut taxon = 0;
        for _ in 0..200 {
            let internal = tree.add_node(Some(current), 0.1, None);
            tree.add_node(Some(current), 0.1, Some(taxon));
            taxon += 1;
            current = internal;
        }
        tree.nodes[current].taxon = Some(taxon);
        tree.root = root;
        let pi = model.stationary();
        let states: Vec<Option<usize>> = (0..=taxon).map(|_| Some(1usize)).collect();
        let ll = log_likelihood(&tree, &model, &TipStates { states: &states }, &pi);
        assert!(ll.is_finite(), "deep-tree likelihood underflowed: {ll}");
    }
}
