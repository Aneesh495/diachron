//! Independent brute-force likelihood by enumeration of internal states.
//!
//! This module exists solely to validate `pruning.rs`. It computes the exact
//! same quantity (the probability of observed tip states marginalized over
//! internal node states) by explicitly iterating over every assignment of
//! states to internal nodes and summing the joint probability. It shares no code
//! with the pruning recursion beyond the tree structure and the model's
//! transition matrices, so agreement between the two is meaningful evidence that
//! the fast algorithm is correct (gate G02).
//!
//! Cost is `n_states ^ (number of internal nodes)`, so this is only used on tiny
//! trees. That is exactly the regime the oracle tests need.

use crate::subst::SubstModel;
use crate::tree::Tree;

/// Enumerate the likelihood of a single character. `tip_states[taxon]` is the
/// observed state or None for missing.
pub fn log_likelihood_enumerate(
    tree: &Tree,
    model: &dyn SubstModel,
    tip_states: &[Option<usize>],
    root_prior: &[f64],
) -> f64 {
    let n_states = model.n_states();
    let n_nodes = tree.nodes.len();

    // Collect internal node indices; these are the free variables we sum over.
    let internal: Vec<usize> = (0..n_nodes).filter(|&i| !tree.is_tip(i)).collect();
    let n_internal = internal.len();

    // Precompute transition matrices per node (branch above it).
    let mut trans = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        trans.push(model.transition(tree.nodes[i].branch_len));
    }

    // Iterate over every assignment of states to internal nodes.
    let total_assignments = (n_states as u64).pow(n_internal as u32);
    let mut acc = 0.0f64;

    // `assign[node] = state` for all nodes; tips are filled per term when the
    // tip is observed, and summed over when missing.
    for code in 0..total_assignments {
        // Decode the internal assignment from a base-n_states integer.
        let mut rem = code;
        let mut state_of = vec![usize::MAX; n_nodes];
        for &node in &internal {
            state_of[node] = (rem % n_states as u64) as usize;
            rem /= n_states as u64;
        }

        // For missing tips we must also sum over their states. Collect them.
        let missing_tips: Vec<usize> = (0..n_nodes)
            .filter(|&i| tree.is_tip(i) && tip_states[tree.nodes[i].taxon.unwrap()].is_none())
            .collect();
        let n_missing = missing_tips.len();
        let missing_terms = (n_states as u64).pow(n_missing as u32);

        for mcode in 0..missing_terms {
            let mut mrem = mcode;
            for &tip in &missing_tips {
                state_of[tip] = (mrem % n_states as u64) as usize;
                mrem /= n_states as u64;
            }
            // Observed tips take their fixed state.
            for i in 0..n_nodes {
                if tree.is_tip(i) {
                    if let Some(s) = tip_states[tree.nodes[i].taxon.unwrap()] {
                        state_of[i] = s;
                    }
                }
            }
            // Joint probability of this full assignment: root prior times the
            // product of edge transition probabilities.
            let mut term = root_prior[state_of[tree.root]];
            for i in 0..n_nodes {
                if let Some(p) = tree.nodes[i].parent {
                    let from = state_of[p];
                    let to = state_of[i];
                    term *= trans[i].at(from, to);
                }
            }
            acc += term;
        }
    }

    if acc <= 0.0 {
        f64::NEG_INFINITY
    } else {
        acc.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pruning::{log_likelihood, TipStates};
    use crate::subst::{BinaryGainLoss, F81};
    use crate::tree::Tree;

    fn balanced_four_tip() -> Tree {
        // root -> {i1 -> {t0,t1}, i2 -> {t2,t3}}
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let i1 = tree.add_node(Some(root), 0.3, None);
        let i2 = tree.add_node(Some(root), 0.5, None);
        tree.add_node(Some(i1), 0.2, Some(0));
        tree.add_node(Some(i1), 0.4, Some(1));
        tree.add_node(Some(i2), 0.1, Some(2));
        tree.add_node(Some(i2), 0.6, Some(3));
        tree.root = root;
        tree
    }

    #[test]
    fn pruning_equals_enumeration_binary() {
        let model = BinaryGainLoss::new(0.4, 0.9);
        let tree = balanced_four_tip();
        let pi = model.stationary();
        let patterns = [
            vec![Some(1), Some(1), Some(0), Some(1)],
            vec![Some(0), Some(0), Some(0), Some(0)],
            vec![Some(1), None, Some(1), Some(0)],
            vec![None, None, Some(1), Some(1)],
        ];
        for pat in &patterns {
            let a = log_likelihood(&tree, &model, &TipStates { states: pat }, &pi);
            let b = log_likelihood_enumerate(&tree, &model, pat, &pi);
            assert!((a - b).abs() < 1e-10, "pattern {pat:?}: {a} vs {b}");
        }
    }

    #[test]
    fn pruning_equals_enumeration_f81() {
        let model = F81::new(vec![0.25, 0.25, 0.5], 1.1);
        let tree = balanced_four_tip();
        let pi = model.stationary();
        let pat = vec![Some(0), Some(2), Some(1), None];
        let a = log_likelihood(&tree, &model, &TipStates { states: &pat }, &pi);
        let b = log_likelihood_enumerate(&tree, &model, &pat, &pi);
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }
}
