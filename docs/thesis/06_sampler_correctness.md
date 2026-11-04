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

Each move proposes a uniform increment. For an internal non-root node,
the window is a fraction of the current root age, which the move leaves
fixed. For the root, the window depends on its oldest child. A dated tip
uses a fixed window from its date density. A proposal that crosses a child
or parent age is rejected. Because the window is unchanged by an accepted
move, the reverse increment has the same density.

## Root scale move

All internal node ages are multiplied by a common factor `exp(lambda(u-1/2))`.
Scaling k independent coordinates contributes a Jacobian of `factor^k`, so the
Hastings term is `k log(factor)`. Invalid results (negative branch) are reverted
and rejected.

## Topology moves

**NNI.** A node's sibling is swapped with one of the node's children. Under
symmetric selection the move is its own inverse, so the Hastings term is zero. The
move is rejected if it would make a child older than its parent.

**Time-preserving SPR.** Choose a prunable child. Detach it and collapse
its parent, keeping the parent's age as a fixed proposal coordinate.
Choose uniformly among edges of the detached backbone crossing that age.
The original sibling edge remains available as a self-transition. In the
reverse direction the same child produces the same detached backbone and
age, so the target-edge sets and selection probabilities agree. The move
preserves the number of prunable children and has no height Jacobian.

The ranked-tree test enumerates all 18 four-tip states at three fixed
event ages. Production NNI and SPR proposals with the exact lexical
target visit all 18 states and meet a fixed frequency tolerance.
This conditional test does not establish mixing at full study scale.

## Checkpointing and exact resumption

The chain's RNG is a serializable xoshiro256** with a four-word state. A
checkpoint writes the RNG words, tree, parameters, frozen step sizes, iteration
counters, and retained draws atomically via a temp-file rename. Resuming
reconstructs an identical chain. A chain interrupted and resumed partway produces
a trajectory bit-identical to an uninterrupted run: same final log posterior,
same topology, identical retained draws, and no duplicated retained draw (gate G12,
`tests/test_recovery.py`). Adaptation of step sizes runs only during warmup and is
frozen before any retained sample, so the retained chain targets the stationary
posterior.
