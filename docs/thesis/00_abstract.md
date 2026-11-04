# Diachron: posterior inference over language family history under uncertain evidence

## Abstract

Diachron is a research platform for studying how language families diversify,
borrow words, and shift their geographic footprint, built around original
probabilistic inference rather than narrative. The core computes a Felsenstein
pruning likelihood of cognate-class evolution under a continuous-time Markov
model, composes it with clock, calibration, tree-prior, geographic, and bounded
borrowing terms into an explicit posterior, and samples that posterior with an
original Metropolis-Hastings kernel whose moves carry derived and tested proposal
ratios and Jacobians.

The thesis makes two kinds of statement and keeps them apart. Mathematical
statements concern the algorithms: the pruning recursion computes the stated
marginal likelihood, and the sampler leaves the stated target invariant. These
are proved against independent tiny-state enumeration (agreement within 1e-9 over
1200 random cases) and against exact finite posteriors. Historical statements are
empirical claims with uncertainty: a posterior over a root age or a latent
homeland is a distribution with credible regions, prior sensitivity, and
counterevidence, never a recovered fact.

A central finding, established on labeled-synthetic known-truth corpora, is that
the inferred chronology is identifiable only through calibration: the root-age
estimate shifts materially when calibrations are removed or halved, while it is
stable to the tree-prior strength and the ascertainment correction. Posterior
predictive checks show the tree model reproduces the pairwise cognate-sharing
structure of the data but does not reproduce the missingness pattern, which the
likelihood handles but the generative check does not, a misspecification we report
rather than hide. Geographic conclusions are presented as posterior clouds; the
platform refuses to collapse a broad or multimodal location posterior to a single
homeland pin, and separates modern reference coordinates, dated attestations, and
inferred latent locations in every view and export.

The contribution is an auditable, reproducible inference pipeline with honest
uncertainty, not a claimed reconstruction of any particular family's history.
