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

/// Conditional probability of the actual binary ascertainment event, given
/// the observed-tip mask. A bucket tracks whether we saw any zero and how
/// many ones we saw (saturated at two). The final accepted bucket is read
/// directly, so near-constant trees do not lose precision to subtraction.
fn retained_probability(
    plan: &PruningPlan<'_>,
    multiplicity: &[usize],
    min_present: usize,
    prior: &[f64],
) -> f64 {
    const BUCKETS: usize = 6;
    const ANY_ZERO: usize = 3;
    type Counts = [f64; BUCKETS];

    let mut messages = vec![[[0.0; BUCKETS]; 2]; plan.tree.nodes.len()];
    for &node in &plan.order {
        let mut mass = [[0.0; BUCKETS]; 2];
        if let Some(taxon) = plan.tree.nodes[node].taxon {
            let copies = multiplicity[taxon];
            mass[0][if copies > 0 { ANY_ZERO } else { 0 }] = 1.0;
            mass[1][copies] = 1.0;
        } else {
            mass[0][0] = 1.0;
            mass[1][0] = 1.0;
            for &child in &plan.tree.nodes[node].children {
                let p = &plan.transitions[child];
                for (state, left) in mass.iter_mut().enumerate() {
                    let mut right: Counts = [0.0; BUCKETS];
                    for (bucket, prob) in right.iter_mut().enumerate() {
                        *prob = p.at(state, 0) * messages[child][0][bucket]
                            + p.at(state, 1) * messages[child][1][bucket];
                    }
                    let mut combined: Counts = [0.0; BUCKETS];
                    for (a, &pa) in left.iter().enumerate() {
                        if pa == 0.0 {
                            continue;
                        }
                        for (b, &pb) in right.iter().enumerate() {
                            let has_zero = (a >= ANY_ZERO || b >= ANY_ZERO) as usize;
                            let ones = (a % ANY_ZERO + b % ANY_ZERO).min(2);
                            combined[has_zero * ANY_ZERO + ones] += pa * pb;
                        }
                    }
                    *left = combined;
                }
            }
        }
        messages[node] = mass;
    }
    messages[plan.tree.root]
        .iter()
        .zip(prior)
        .map(|(mass, pi)| pi * mass[ANY_ZERO + min_present..].iter().sum::<f64>())
        .sum()
}

fn log_mean_exp(log_probs: &[f64]) -> f64 {
    let max = log_probs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return max;
    }
    let total: f64 = log_probs.iter().map(|p| (p - max).exp()).sum();
    max + (total / log_probs.len() as f64).ln()
}

/// Equal-weight site rates and a finite contact mixture with baseline weight 1.
/// At most one edge copies per character. Ascertainment conditions the mixture
/// after summing both the pattern mass and selection mass over components.
pub fn binary_mixture_scores(
    tree: &Tree,
    model: &dyn SubstModel,
    columns: &[Vec<Option<usize>>],
    root_prior: &[f64],
    rate_categories: &[f64],
    min_present: usize,
    contacts: &[(usize, usize, f64)],
) -> Vec<f64> {
    use std::collections::HashMap;

    let plans: Vec<_> = rate_categories
        .iter()
        .map(|rate| PruningPlan::new(tree, model, *rate))
        .collect();
    let mut normalizers = HashMap::new();
    let mut mixtures = vec![0.0; plans.len()];
    let weight_sum = 1.0 + contacts.iter().map(|(_, _, w)| w).sum::<f64>();
    let mut working_column = vec![None; tree.n_tips()];
    let mut components = vec![0.0; contacts.len() + 1];
    let score = |column: &[Option<usize>], mixtures: &mut [f64]| {
        for (result, plan) in mixtures.iter_mut().zip(&plans) {
            *result = plan.score(column, root_prior);
        }
        log_mean_exp(mixtures)
    };
    columns
        .iter()
        .map(|column| {
            let log_normalizer = if min_present == 0 {
                weight_sum.ln()
            } else {
                let present = column.iter().filter(|s| **s == Some(1)).count();
                if present < min_present || !column.contains(&Some(0)) {
                    return f64::NEG_INFINITY;
                }
                let observed: Vec<_> = column
                    .iter()
                    .map(|state| usize::from(state.is_some()))
                    .collect();
                *normalizers.entry(observed.clone()).or_insert_with(|| {
                    let retained = |multiplicity: &[usize]| {
                        plans
                            .iter()
                            .map(|plan| {
                                retained_probability(plan, multiplicity, min_present, root_prior)
                            })
                            .sum::<f64>()
                    };
                    let base = retained(&observed);
                    let mut sum = base;
                    let mut counts = observed.clone();
                    for &(donor, recipient, intensity) in contacts {
                        if observed[recipient] == 0 {
                            sum += intensity * base;
                        } else {
                            counts[donor] += 1;
                            counts[recipient] = 0;
                            sum += intensity * retained(&counts);
                            counts[donor] -= 1;
                            counts[recipient] = 1;
                        }
                    }
                    (sum / plans.len() as f64).ln()
                })
            };
            if !log_normalizer.is_finite() {
                return f64::NEG_INFINITY;
            }
            let base = score(column, &mut mixtures);
            components[0] = base;
            for (&(donor, recipient, intensity), component) in
                contacts.iter().zip(&mut components[1..])
            {
                let log_mass = match (column[donor], column[recipient]) {
                    (_, None) => base,
                    (Some(d), Some(r)) if d != r => f64::NEG_INFINITY,
                    (_, observed) => {
                        working_column.copy_from_slice(column);
                        working_column[donor] = observed;
                        working_column[recipient] = None;
                        score(&working_column, &mut mixtures)
                    }
                };
                *component = log_mass + intensity.ln();
            }
            let log_weighted_sum = log_mean_exp(&components) + (components.len() as f64).ln();
            log_weighted_sum - log_normalizer
        })
        .collect()
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
