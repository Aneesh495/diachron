# Note 5: Clock and date models

A clock turns years into expected substitutions, and a calibration anchors the
years to the calendar. The two together are what let a cognate matrix speak to
chronology, and understanding how they interact is the key to reading the results
honestly.

The strict clock has one rate for the whole tree. Without a calibration, the rate
and the absolute time scale trade off exactly: doubling the rate and halving all
ages leaves the substitution counts, and hence the likelihood, unchanged. This is
not a bug; it is the identifiability structure of the problem. A calibration
breaks the symmetry by pinning a node age to an interval, after which the rate and
the ages are jointly identified.

Diachron derives one aggregate root-age calibration from the dated attestations
rather than one calibration per attestation. A single contemporary tip at age zero
carries no time-scale information, and placing a tight age constraint on every
tip's parent over-constrains the tree to the point that nothing is accepted (an
early version did exactly this and the sampler stalled at zero acceptance). The
aggregate root calibration is the minimal structure that identifies the time scale
while leaving the topology free to move.

The relaxed uncorrelated-lognormal clock gives each branch its own rate and
reduces to the strict clock as its log-standard-deviation goes to zero. Whether
the relaxed clock is identifiable from a given dataset is itself a question: with
sparse calibration the branch-rate multipliers are weakly constrained, and the
honest report is a posterior that is close to the prior, not a confident claim of
rate heterogeneity.
