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

Products of many small probabilities underflow. At each internal node the
conditional vector is divided by its maximum entry, and the log of that maximum is
accumulated into a per-node log scale; the final log likelihood adds the
accumulated scales to `log sum_s pi(s) L_root(s)`. This keeps the stored vector
O(1) regardless of tree depth. A caterpillar tree of 200 tips, which underflows
without scaling, returns a finite log likelihood with it
(`tests/` and `cli stress`).

## Missing versus absent

A missing tip marginalizes over its state; an absent tip pins the state to 0.
Because marginalizing sums the two indicator cases, the missing likelihood is
strictly larger than either pinned case. The Rust and Python tests both confirm
`L(missing) > L(absent)` for a fixed surrounding pattern (gate G04). This is the
mechanism by which "we have no data here" and "the class is absent here" produce
different likelihoods, as the evidence contract requires.

## Ascertainment conditioning

Only variable characters are kept. Conditioning each kept character on being
variable divides its likelihood by `p(variable) = 1 - p(constant)`, where
`p(constant)` is the summed probability of the all-absent and all-present patterns
under the same tree and model. We compute `p(constant)` by running the pruning
likelihood on the two constant patterns and summing. Validation: the analytic
`p(variable)` matches the empirical fraction of forward simulations that come out
variable, to within 0.005 over 400000 simulations (gate G04,
`tests/test_ascertainment.py`). The correction matches the selection process
actually applied; a correction for a different scheme would not reproduce the
simulated fraction.

## Geographic Gaussian pruning

The Brownian-on-tree likelihood uses independent contrasts: at each internal node
the children's location estimates are combined by precision weighting, the added
variance is the harmonic combination of the child branch variances, and each
child-pair contrast contributes a Gaussian density. This integrates out all
internal locations in closed form and gives the exact marginal likelihood of the
tip coordinates under Brownian motion, the continuous counterpart of the discrete
pruning recursion.
