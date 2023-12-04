# Note 6: Proposal correctness

A Metropolis-Hastings move is only correct if its Hastings ratio is right, and the
two places that is easy to get wrong are the Jacobian of a transformed-variable
move and the proposal ratio of a dimension- or structure-changing move.

The scale move multiplies a positive scalar by `exp(lambda(u-1/2))`. Its proposal
is symmetric in log space, so the Hastings term is the Jacobian of the transform,
`log(factor)`. The test that this is right is not that a Gamma target is recovered
with the Jacobian; it is that a negative control with the Jacobian omitted
recovers a detectably wrong mean. Without the control, a passing Gamma test could
be vacuous (the move might be nearly symmetric enough to pass by luck); with it,
the Jacobian is shown to be doing real work.

The root-scale move multiplies k internal node ages by a common factor, so its
Jacobian is `factor^k` and the Hastings term is `k log(factor)`. The topology
moves are the subtler case. NNI under symmetric selection is its own inverse, so
its Hastings term is zero. Bounded SPR is not symmetric: the number of valid
regraft targets differs before and after, so the Hastings ratio is the log of
that count ratio. Getting the valid-target set right is also a correctness matter
beyond the ratio: an earlier SPR omitted the exclusions that prevent reattaching
inside the pruned subtree, which created cycles. The cycle did not just bias the
sampler; it made the likelihood traversal loop, turning a two-second chain into a
ninety-second one, which is how the bug was found. The lesson recorded here is
that a structure-changing move needs an integrity test (no cycles, one root,
stable node count) in addition to its ratio derivation.

## Coupled rate and clock proposals

For log step d symmetric about zero, transform
(g, l, c) to (g exp(-d), l exp(-d), c exp(d)). The two effective
substitution rates g c and l c are invariant. In original coordinates the
auxiliary-variable transform has absolute Jacobian exp(-d). The log
Hastings term is -d. An independent test conditions three lognormal priors
on fixed g c and l c, derives the analytic normal law for log(c), and
checks the production proposal against it.

The joint age and clock move multiplies k internal ages by exp(d) and
divides c by exp(d), with log Jacobian (k-1)d. Contemporary tip branches
preserve their substitution lengths. Historical tip branches change, and
the full posterior evaluates those changes. Invalid time order is rejected.
