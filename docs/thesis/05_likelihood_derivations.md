# Chapter 5: Likelihood derivations and numerical stability

## Pruning recursion correctness

Claim (mathematical). For a rooted tree with branch lengths, a reversible
two-state CTMC with transition matrices P(t), and observed tip states, the pruning
recursion computes exactly

    P(observed tips) = sum over all internal-state assignments of
                       pi(root state) * prod_edges P(parent state -> child state).

Proof sketch. The recursion is the distributive law applied to the sum over
internal states. Factor the joint probability as a product over edges; each
internal node's state appears only in the edges incident to it, so the sum over
that state can be pushed inward past all terms not involving it. Processing nodes
in post-order performs exactly these inner sums, yielding `L_v(s)` as the partial
sum over the subtree below v conditional on v's state. At the root the remaining
sum against pi is the full marginal. A missing tip contributes the all-ones vector
because its unobserved state is summed over all values, which is the marginal over
that tip; an observed tip contributes an indicator.

Verification (empirical). An independent enumerator
(`rust/diachron-core/src/enumerate.rs`, and a separate Python oracle in
`tests/oracle.py`) sums the joint probability over all internal-state assignments
with no shared code path. Over 1200 random trees with 2 to 5 tips, random rates,
and 15% missing tips, the pruning and enumeration log likelihoods agree to within
1e-9 (gate G02, `tests/test_likelihood_oracle.py`). The agreement of two
independently written routines in two languages is the evidence that the fast
algorithm is correct.

## Log-space scaling

Products of many small probabilities underflow. After each child contribution the
conditional vector is divided by its maximum entry, and the log of that maximum is
accumulated into a per-node log scale; the final log likelihood adds the
accumulated scales to `log sum_s pi(s) L_root(s)`. This keeps the stored vector
O(1) regardless of tree depth. A caterpillar tree of 200 tips, which underflows
without scaling, returns a finite log likelihood with it
(`tests/` and `cli stress`).

## Missing versus absent

A missing tip marginalizes over its state; an absent tip pins the state to 0.
Because marginalizing sums the two indicator cases, the missing likelihood is
at least as large as either pinned case. Equality is possible if the other state has zero mass. The Rust and Python tests both confirm
`L(missing) > L(absent)` for a fixed surrounding pattern (gate G04). This is the
mechanism by which "we have no data here" and "the class is absent here" produce
different likelihoods, as the evidence contract requires.

## Ascertainment conditioning

Condition on retention given the observed mask. The default rule needs
at least two observed presences and at least one observed absence.
When singletons are allowed, the positive threshold is one. Unobserved
tips contribute zero to both counts.

For each ancestral state, the normalization pass tracks six probabilities.
A bucket records the presence count saturated at two and a flag for any
absence. The transition step propagates child buckets through the CTMC.
The merge convolves saturated counts and takes the logical union of
absence flags. At the root the stationary-weighted mass of acceptable
buckets is the retention probability. For an observed recipient in a
contact copy, count the donor twice if both tips are observed, and once
if only the recipient is observed. The recipient's inherited state then
contributes no count. Rate and contact probabilities are averaged before
taking the ratio. The tests enumerate retained patterns with both
thresholds, varied masks, and contact components.

## Geographic Gaussian pruning

The Brownian-on-tree likelihood uses independent contrasts: at each internal node
the children's location estimates are combined by precision weighting, the added
variance is the harmonic combination of the child branch variances, and each
child-pair contrast contributes a Gaussian density. This integrates out all
internal locations in closed form and gives the exact likelihood of contrasts under the stated planar
Gaussian model, the continuous counterpart of the discrete
pruning recursion.
