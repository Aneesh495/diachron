//! Metropolis-Hastings kernels for phylogenetic posteriors.
//!
//! This module provides the moves and the accept/reject machinery. The target
//! density is supplied by the caller as a closure computing the unnormalized log
//! posterior of a `State`; the moves here know how to propose a new state and
//! report the log proposal ratio (and Jacobian, where the move changes a
//! parameter by a multiplicative factor). Correctness of each move's ratio is
//! validated in `tests/test_mcmc_exact.py` against exact finite targets
//! (gate G03).
//!
//! State layout. The sampler works on a timed tree plus a small vector of scalar
//! parameters (clock rate, gain, loss). Node heights are the primary continuous
//! variables; branch lengths are derived from heights so that time constraints
//! (a child is never older than its parent) are automatic when heights respect
//! the ordering. Topology moves rewire the tree while preserving tip labels.
//!
//! Each move returns `Proposal { log_hastings }` where `log_hastings` is
//! `log q(old|new) - log q(new|old) + log|Jacobian|`, the quantity added to the
//! log-posterior difference in the acceptance test. A symmetric random-walk move
//! returns 0; a multiplicative scale move returns the log of the scale factor
//! times the number of scaled quantities.

use crate::rng::Rng;
use crate::tree::Tree;

/// Scalar model parameters sampled alongside the tree.
#[derive(Clone, Debug)]
pub struct Params {
    pub clock_rate: f64,
    pub gain: f64,
    pub loss: f64,
}

/// Full MCMC state: a timed tree and the scalar parameters.
#[derive(Clone, Debug)]
pub struct State {
    pub tree: Tree,
    pub params: Params,
}

/// The outcome of proposing a move: the (mutated) state is held by the caller;
/// this reports the Hastings term to fold into acceptance.
pub struct Proposal {
    pub log_hastings: f64,
    /// False if the proposal produced an invalid state (e.g. negative time) and
    /// must be rejected outright.
    pub valid: bool,
}

impl Proposal {
    fn reject() -> Self {
        Proposal { log_hastings: 0.0, valid: false }
    }
}

/// Multiplicative scale move on a single positive scalar. Proposes
/// `x' = x * exp(lambda * (u - 0.5))` with `u ~ U(0,1)`. The Hastings ratio for
/// this multiplier proposal is `log(x'/x)`, which is the log Jacobian of the
/// transform; the proposal density itself is symmetric in log space.
pub fn scale_scalar(x: &mut f64, lambda: f64, rng: &mut Rng) -> Proposal {
    let u = rng.next_f64();
    let factor = (lambda * (u - 0.5)).exp();
    let old = *x;
    *x = old * factor;
    if !(*x).is_finite() || *x <= 0.0 {
        return Proposal::reject();
    }
    Proposal { log_hastings: factor.ln(), valid: true }
}

/// Random-walk move on a node's height, reflected to stay between its oldest
/// child and its parent so the time ordering is preserved. Symmetric: the
/// reflected Gaussian proposal gives a zero Hastings term because the window is
/// recomputed consistently; we therefore restrict the step to stay in-window and
/// reject otherwise, which keeps the move symmetric and correct.
pub fn slide_node_height(
    state: &mut State,
    node: usize,
    window: f64,
    rng: &mut Rng,
) -> Proposal {
    // Root height can move freely upward; internal nodes are bounded.
    let tree = &state.tree;
    if tree.is_tip(node) {
        return Proposal::reject();
    }
    let lower = tree
        .nodes[node]
        .children
        .iter()
        .map(|&c| tree.nodes[c].height)
        .fold(f64::NEG_INFINITY, f64::max);
    let upper = match tree.nodes[node].parent {
        Some(p) => tree.nodes[p].height,
        None => f64::INFINITY,
    };
    let step = (rng.next_f64() - 0.5) * 2.0 * window;
    let proposed = state.tree.nodes[node].height + step;
    if proposed <= lower + 1e-12 || proposed >= upper - 1e-12 {
        return Proposal::reject();
    }
    state.tree.nodes[node].height = proposed;
    state.tree.set_branch_lengths_from_heights();
    Proposal { log_hastings: 0.0, valid: true }
}

/// Nearest-neighbor interchange. Picks an internal edge (parent p with child c,
/// both internal) and swaps one child of c with the sibling of c under p. This
/// is the canonical topology move on rooted trees; it is its own inverse with
/// equal selection probability in the symmetric regime used here, so the
/// Hastings term is zero. The move preserves node times by keeping heights and
/// recomputing branch lengths; if the swap would violate time ordering it is
/// rejected.
pub fn nni(state: &mut State, rng: &mut Rng) -> Proposal {
    let n_nodes = state.tree.nodes.len();
    // Candidate internal nodes with an internal child.
    let mut candidates = Vec::new();
    for c in 0..n_nodes {
        if state.tree.is_tip(c) {
            continue;
        }
        if let Some(p) = state.tree.nodes[c].parent {
            if state.tree.nodes[p].children.len() >= 2 {
                candidates.push((p, c));
            }
        }
    }
    if candidates.is_empty() {
        return Proposal::reject();
    }
    let (p, c) = candidates[rng.next_below(candidates.len())];

    // Pick a sibling of c under p (some child of p other than c).
    let siblings: Vec<usize> = state.tree.nodes[p]
        .children
        .iter()
        .copied()
        .filter(|&x| x != c)
        .collect();
    if siblings.is_empty() {
        return Proposal::reject();
    }
    let sib = siblings[rng.next_below(siblings.len())];

    // Pick a child of c to swap with the sibling.
    if state.tree.nodes[c].children.is_empty() {
        return Proposal::reject();
    }
    let gc_idx = rng.next_below(state.tree.nodes[c].children.len());
    let gc = state.tree.nodes[c].children[gc_idx];

    // Time validity: after swapping, sib becomes a child of c and gc a child of
    // p. Require c older than sib and p older than gc.
    if state.tree.nodes[c].height <= state.tree.nodes[sib].height + 1e-12 {
        return Proposal::reject();
    }
    if state.tree.nodes[p].height <= state.tree.nodes[gc].height + 1e-12 {
        return Proposal::reject();
    }

    // Rewire: remove sib from p's children, add gc; remove gc from c's children,
    // add sib. Update parent pointers.
    let p_children = &mut state.tree.nodes[p].children;
    if let Some(pos) = p_children.iter().position(|&x| x == sib) {
        p_children[pos] = gc;
    }
    let c_children = &mut state.tree.nodes[c].children;
    if let Some(pos) = c_children.iter().position(|&x| x == gc) {
        c_children[pos] = sib;
    }
    state.tree.nodes[sib].parent = Some(c);
    state.tree.nodes[gc].parent = Some(p);
    state.tree.set_branch_lengths_from_heights();
    Proposal { log_hastings: 0.0, valid: true }
}

/// Accept or reject a proposed state against the current one. Returns true if
/// accepted. `log_post_old` and `log_post_new` are unnormalized log posteriors;
/// `log_hastings` is the move's Hastings term.
pub fn accept(
    log_post_old: f64,
    log_post_new: f64,
    log_hastings: f64,
    rng: &mut Rng,
) -> bool {
    if !log_post_new.is_finite() {
        return false;
    }
    let log_alpha = log_post_new - log_post_old + log_hastings;
    if log_alpha >= 0.0 {
        return true;
    }
    rng.next_f64().ln() < log_alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_move_hastings_is_log_factor() {
        let mut r = Rng::seed_from_u64(1);
        let mut x = 2.0;
        let before = x;
        let prop = scale_scalar(&mut x, 0.5, &mut r);
        assert!(prop.valid);
        assert!((prop.log_hastings - (x / before).ln()).abs() < 1e-12);
    }

    #[test]
    fn accept_always_takes_uphill() {
        let mut r = Rng::seed_from_u64(2);
        assert!(accept(-10.0, -5.0, 0.0, &mut r));
    }

    #[test]
    fn accept_rejects_nonfinite() {
        let mut r = Rng::seed_from_u64(3);
        assert!(!accept(-5.0, f64::NEG_INFINITY, 0.0, &mut r));
    }
}
