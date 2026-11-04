# Note 8: Synthetic recovery and coverage

A recovery study answers a question that real data cannot: given that we know the
true tree, how well does inference find it, and is the uncertainty honest?

The design generates datasets under the inference model, varying tree size,
concept count, root age, borrowing, and missingness, and measures clade recovery
(fraction of true clades with posterior support above one half) and credible-region
coverage (fraction of datasets whose true root age lies in the 90 percent
posterior interval). Coverage is the honesty measure: calibrated intervals contain
the truth about 90 percent of the time.

The pilot of five datasets with short chains gives mean clade recovery 0.41 and
90 percent coverage 0.8. The moderate clade recovery reflects the short chains and
the borrowing perturbation; the coverage near target is early evidence the
intervals are not badly over-confident. The full 200-dataset study runs with the
frozen budget.

The caveat recorded here is that recovery under the model's own assumptions is
necessary but not sufficient. A posterior that is narrower than it should be can
still cover the truth on self-generated data while being over-confident under
misspecification. That is why the recovery regime includes borrowing and
missingness the fitted model does not match, and why Note 10 and the posterior
predictive checks probe fit to features the generative model does not target.
