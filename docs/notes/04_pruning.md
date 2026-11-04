# Note 4: The pruning likelihood

The pruning recursion is the one piece of the engine that must be both fast and
exactly right, because everything downstream multiplies through it. Two decisions
make it trustworthy.

First, validation against an independent oracle. The enumerator sums the joint
probability over every internal-state assignment with no shared code path to the
recursion, and the two agree to within 1e-9 over 1200 random tiny trees. Writing
the oracle in a second language (the Rust enumerator and a separate Python one)
guards against a shared transcription bug: an error reproduced in both would have
to be independently reinvented, which an enumeration over states will not do.

Second, log-space scaling. Without rescaling, the conditional likelihood vectors
underflow on deep trees and the log likelihood becomes negative infinity. At each
internal node the vector is divided by its maximum and the log of that maximum is
accumulated; a 200-tip caterpillar that underflows without this returns a finite
value with it. The scaling is exact, not approximate: it factors a positive
constant out of a product and adds its log back, so no precision is traded for
stability.

A measured aside: on a 60-tip tree the binding evaluates roughly ten thousand
character likelihoods per second, which is what makes the Python-driven MCMC
practical despite crossing the language boundary on every step.
