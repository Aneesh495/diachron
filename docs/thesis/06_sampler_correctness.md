# Chapter 6: Sampler correctness

The sampler is a Metropolis-within-Gibbs cycle over scalar parameters and tree
structure. This chapter states the proposal ratios and describes the tests
that compare production moves with exact targets.

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
the window is a fraction of the interval from its oldest child to its
parent. Both bounds remain fixed during that move. For the root, the window depends on its oldest child. A dated tip
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

## Coupled rate and time proposals

The rate-clock ridge proposal multiplies clock by f and divides gain and
loss by f, preserving effective lexical rates. Its log Hastings term is
-log(f), from the product of the original coordinates. The independent
conditional-lognormal test checks the induced clock distribution.

The coupled time-clock move multiplies k internal ages by f and divides
clock by f. Its log Hastings term is (k-1) log(f). Dated tips stay fixed,
so effective branch lengths need not be invariant. The acceptance ratio
reevaluates the complete dated target. The dedicated root move scales the
positive gap above the root's oldest child and contributes log(f).

## Composed finite posterior check

`make sampler-check` enumerates 72 conditional states for each required
observation model, combining ranked four-tip topologies, rate settings and
two historical-tip ages. The reference sums complete binary node states
and contact histories, applies the observed-mask retention event, integrates
a dense Gaussian geographic precision matrix and uses independent SciPy
prior and terminal-date densities. It shares no pruning or geographic
message implementation with the production score.

The production acceptance method samples each finite target for 200000
retained steps after 1000 warmup steps. The predeclared maximum frequency
deviation is 0.01 and the log-posterior tolerance is 1e-9. Passing this audit
validates the composed conditional targets. Continuous proposal ratios,
checkpoint recovery and full-data mixing require their separate checks.

## Resampled attachment age

The optional native SPR move prunes the same kind of child as fixed-height
SPR. It suppresses the binary parent. Each backbone edge with enough time
above the detached child is eligible. The proposal chooses an edge uniformly
and samples the reattached parent's age uniformly between the oldest child
and the edge's upper endpoint. The original edge remains eligible.

Detaching the same child after the move yields the same backbone. The number
of prunable children and the eligible edge set agree in reverse. If the old
and new age spans are w and w', the age density ratio is w'/w. The accepted
proposal contributes log(w') minus log(w). A four-tip test fixes the root
and samples uniformly over the valid internal heights. Each balanced topology
has a square of height volume. Each pectinate topology has half a square.
Sampled frequencies meet the declared tolerance for all fifteen labeled
rooted topologies.

## Weighted attachment candidates

An optional attachment move draws four candidates on the detached backbone
and adds the current state to that set. Given the detached child and
backbone, it selects an entry with weight proportional to its full posterior
density divided by its edge and age proposal density. In the augmented
distribution, the selected position has the target density and the other
positions have the proposal density. Resampling the position from those
weights preserves the target marginal.

The parsimony option uses two-state inside and outside messages to compute
the minimum lexical changes at each potential insertion edge. Those counts
use only the fixed columns and the shared backbone and subtree. Missing
states allow either value. The finite-edge probabilities retain a uniform
component and the root stem retains its own mass. Enumerating states on
small trees checks the placement scores. A separate four-tip test checks
the stationary topology frequencies and unbounded root-age distribution
with and without the parsimony weights.

## Joint calendar scale and time block

The calendar-scale move multiplies the inferred ages by exp(t) and divides
the clock and any active diffusion rate by exp(t). Its volume factor is
exp((d-r)t) for d free ages and r reciprocally scaled rates. Random placement
of a unit interval on the log scale and subsequent shrinkage
produce a reversible slice kernel for the calendar, coalescent and rate
priors times that volume factor. The full posterior acceptance calculation
subtracts the prior-density change accounted for by the slice kernel.
[Slice sampling](https://arxiv.org/abs/physics/0009028)

The joint time and clock block uses log gaps above the oldest child for
internal nodes and log gaps above the date offset for sampled tips. The log
clock is a separate coordinate. Computing ages from younger to older gives
a triangular change of variables with log determinant equal to the sum of
those log coordinates. A Gaussian proposal uses the lexical, coalescent and
terminal date gradients to guide its mean. Acceptance uses the full posterior,
both Gaussian densities and the log determinant difference. The step size
adapts during warmup and is saved in the checkpoint.
[Langevin proposal reference](https://doi.org/10.2307/3318418)

With multiple steps the same chart supports a Hamiltonian proposal. Fixed
metric drifts and position-dependent momentum kicks each preserve volume.
The symmetric sequence followed by a momentum flip is its own inverse.
The acceptance rule adds the initial minus final Gaussian kinetic energy to
the full posterior and chart determinant changes. This remains valid when
the deterministic force omits a posterior term. A test checks the inverse
even for a nonconservative force, and a separate run resumes the complete
multi-step kernel exactly.

The lexical gradient combines derivatives of the observed pattern and the
retention probability. The pattern derivative pairs a descendant likelihood
with its complement on each edge. The retention derivative carries the
same saturated count states as the likelihood, along with outside weights
for each count bucket. Rate and contact mixtures average component
derivatives by their mass. Tests compare the result with finite differences
of the production likelihood and check the time-coordinate determinant
and exact resume. The observed pilots measure mixing separately.
[Sampler pilot results](../results/sampler-pilots.md)

## Lexical replica exchange

For each replica with inverse temperature beta, its target is proportional
to the lexical likelihood to the power beta times the other posterior factors.
The temperature-one marginal is the original posterior. The joint target is
the product across replicas. Swapping two complete states leaves the product
of all nonlexical factors fixed. With lexical scores L_a and L_b, the log
exchange ratio is (beta_a - beta_b)(L_b - L_a). Each local transition and each
adjacent exchange preserves the joint target after warmup.

A finite test computes the entire two-replica joint distribution directly and
compares its probabilities with sampled frequencies. The production checkpoint
stores all local RNGs, the exchange RNG, the current walker order, swap counts
and cold draws in one atomic file. An interrupted run resumes exactly at a
nonexchange iteration, including all hot states and cold retained iterations.
