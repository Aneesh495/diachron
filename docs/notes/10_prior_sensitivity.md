# Note 10: Prior sensitivity

Fifteen configurations completed the protocol 2 synthetic diagnostic at
`11fe5e485071f3984dec09090525354d89a0e4da`. Each used two chains, 500 warmup
iterations and 1500 sampling iterations. The exact configurations, root
samples' quantiles and descriptive shifts are in
[the saved report](../results/sensitivity-v2.json).

This diagnostic does not establish robustness. In particular, the
coalescent-size settings produced substantially different root means. Short
chains, distinct topology modes and changes to the prior can all cause those
shifts. The experiment does not isolate their contributions. An earlier
version of this note claimed that the estimate was stable across these tree
priors. The saved result does not support that claim.

The calibration subset removes alternate synthetic attestation records before
deriving a root constraint. It can preserve the oldest attestation, so it does
not guarantee a weaker root constraint. For observed IE-CoR, dating evidence
is attached to terminal varieties. EDH inscriptions have no resolved mapping
to these tips and do not calibrate the posterior. Halving the EDH list would
not test sensitivity to IE-CoR dates.

An observed calibration intervention must state which terminal densities
change and which replacement prior keeps a historical tip's age latent. Loan
removal must propagate uncertain absences to missing observations. Alternative
cognacy assignments require an explicit transformation of source judgments.
Every intervention needs its own saved model and source reconstruction.

The decision rule remains: a historical conclusion is stable only if relevant
fits converge and their changes remain within the declared comparison scope.
Small differences between unconverged fits are insufficient evidence of
robustness. Large differences are a reason to investigate the prior, data
transformation and mixing separately.
