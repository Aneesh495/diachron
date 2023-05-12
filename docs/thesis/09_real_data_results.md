# Chapter 9: Results on the study corpus

Results here are on the primary study corpus (112 varieties, 110 concepts, 1576
variable characters, 140 attestations). The corpus is labeled synthetic
known-truth; the same pipeline runs on CLDF collections ingested under their own
licenses.

## Inference and convergence

A diagnostic four-chain run of the dated model (1200 retained iterations per
chain after 400 warmup) produces a root-age posterior with mean about 7240 years
and active topology exploration: 60 to 72 distinct topologies per chain, a
topology switch rate near 0.55, and a between-chain clade agreement of 0.76. The
scalar endpoints on this short run have split-Rhat between 1.1 and 2.8 and ESS in
the tens, so the run is reported as not converged; it is a diagnostic profile, not
the acceptance campaign. The honest reading is that a short run of a 112-tip tree
space has not mixed, and more iterations are required before the posterior is
reported as converged.

## Held-out prediction

Held-out-concept prediction holds out whole concepts so no cognate class spans
the partition, fits on the training characters only, and scores the held-out
characters by posterior predictive log density per character. The comparison
between the dated and tree-only models, and the leakage guard, are recorded in
`runs/predict/heldout.json`. The design, not a single number, answers Q1:
geographic or dated structure earns its place only if it improves held-out
density over the tree-only baseline.

## Posterior predictive checks

Simulating characters from posterior draws and comparing summary statistics to
the observed data gives a Bayesian p-value near 0.18 for mean pairwise cognate
sharing, so the tree model reproduces who-shares-cognates well. The missingness
rate has a p-value of 0 because the generative check does not model the
missingness process that the likelihood handles by marginalization; this is a
reported misspecification of the generative check, not of the likelihood. The
mean class size is slightly over-predicted (p near 0.93). These results are the
kind of honest mixed fit the platform is built to surface.
