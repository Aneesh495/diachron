# Chapter 6: Sampler correctness

The sampler is a Metropolis-within-Gibbs cycle over scalar parameters and tree
structure. This chapter states each move's proposal ratio, derives the two that
carry a Jacobian, and reports the tests that validate them against exact targets.

## Acceptance rule

A proposed state is accepted with probability `min(1, exp(d))` where

    d = log p(proposed) - log p(current) + log_hastings,

and `log_hastings = log q(current|proposed) - log q(proposed|current) + log|J|`.
Non-finite proposals are rejected. The rule is validated on a finite five-state
target: sampled frequencies match the exact normalized probabilities to a maximum
absolute deviation below 0.01 over 400000 iterations (gate G03,
`tests/test_mcmc_exact.py`).

## Scalar scale move

For a positive scalar the move proposes `x' = x * exp(lambda (u - 1/2))`,
`u ~ U(0,1)`. The proposal density is symmetric in log space, so the Hastings
ratio is the Jacobian of the transform, `log(x'/x) = log(factor)`. This is tested
two ways. First, with the correct Jacobian the move targets a Gamma(3, 2) density
with sampled mean 6 and variance 12 (analytic values) within predeclared
tolerances. Second, a negative control runs the same sampler with the Jacobian
omitted and confirms the stationary mean is detectably wrong, so the passing test
is not vacuous.

## Node-time slide

An internal node's age is moved by a uniform step within a window, reflected to
stay strictly between its oldest child and its parent. Because the window is
recomputed consistently and out-of-window proposals are rejected, the move is
symmetric and the Hastings term is zero. Branch lengths are derived from node
ages, so time ordering is maintained automatically.

## Root scale move

All internal node ages are multiplied by a common factor `exp(lambda(u-1/2))`.
Scaling k independent coordinates contributes a Jacobian of `factor^k`, so the
Hastings term is `k log(factor)`. Invalid results (negative branch) are reverted
and rejected.

## Topology moves

**NNI.** A node's sibling is swapped with one of the node's children. Under
symmetric selection the move is its own inverse, so the Hastings term is zero. The
move is rejected if it would make a child older than its parent.

**Bounded SPR.** A node and its parent edge are pruned; the collapsed parent is
reinserted on a target edge whose time interval spans the pruned node's crown age.
Targets are restricted to edges outside the pruned subtree and excluding the
collapsed node and the sibling edge, which is what prevents cycles; the Hastings
ratio is the ratio of valid-target counts before and after. An earlier version
omitted these exclusions and could create a cycle, which made likelihood traversal
blow up; the fix and a 3000-move integrity test that checks for cycles, a single
root, a stable node count, and valid times are in `tests/test_moves.py`.

## Checkpointing and exact resumption

The chain's RNG is a serializable xoshiro256** with a four-word state. A
checkpoint writes the RNG words, tree, parameters, frozen step sizes, iteration
counters, and retained draws atomically via a temp-file rename. Resuming
reconstructs an identical chain. A chain interrupted and resumed partway produces
a trajectory bit-identical to an uninterrupted run: same final log posterior to
1e-9, same topology, same draw count, and no duplicated retained draw (gate G12,
`tests/test_recovery.py`). Adaptation of step sizes runs only during warmup and is
frozen before any retained sample, so the retained chain targets the stationary
posterior.
