# Chapter 1: Research questions and hypotheses

These questions and hypotheses were frozen before the full experiments were run
(see `docs/EXPERIMENT_PLAN.md`). They separate what the evidence can identify
from what it cannot, and each resolves to a concrete comparison in a later
chapter.

## Predeclared questions

**Q1. Geographic versus purely phylogenetic predictive value.** Does conditioning
on geography improve held-out predictive density over a tree-only model, or does
the tree alone carry the signal? Identifiability concern: geographic tip
coordinates are modern reference points, so any geographic gain may reflect
present-day areal structure rather than history.

**Q2. Borrowing sensitivity.** When a known loan set is present, does a
borrowing-aware model recover the contact signal, and does removing loans from the
inherited characters change the tree? Limit: the borrowing model uses a bounded
candidate network, so it cannot recover unrestricted reticulation.

**Q3. Ascertainment bias.** Does the variable-character ascertainment correction
change inferences materially, and does switching missing-to-absent produce the
distinct likelihoods the model requires?

**Q4. Calibration dependence of chronology.** How much does the inferred root age
depend on the calibration priors versus the character likelihood? Hypothesis H4:
without calibration the time scale is unidentifiable and the root age is driven by
the prior.

**Q5. Rate heterogeneity.** Does a relaxed clock change divergence-time estimates
relative to a strict clock, and is the relaxed clock identifiable from this data?

**Q6. Prior sensitivity of geography.** Does the geographic prior (spherical
diffusion variance, or region granularity) drive the latent-location posterior
more than the data?

**Q7. Posterior multimodality.** Do independent chains agree on the tree
posterior, or is there multimodality that scalar diagnostics would miss?

**Q8. Influence of single observations.** Does removing one influential language
or one attestation move an apparently stable homeland or age estimate?

## Strongest initial hypothesis and its counterargument

The strongest initial hypothesis a reader might hold is that a coherent, tightly
dated homeland can be read off the posterior. Chapter 13 argues against this
directly: the platform's own sensitivity results show the chronology is
calibration-dependent and the location posterior is often broad, so a single
confident homeland claim is not supported by this evidence even when a point
estimate looks plausible. Whether that counterargument ultimately wins for a given
family is itself an empirical question the platform is built to ask, not assume.
