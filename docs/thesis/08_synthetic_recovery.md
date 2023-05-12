# Chapter 8: Synthetic recovery

Known-truth recovery answers two questions: can the sampler recover a tree it
generated under its own model, and are its credible regions honest?

## Design

Each synthetic dataset is drawn with a coalescent tree, cognate-class evolution
under the inference model, optional borrowing, Brownian geography, and dated
attestations, with the true tree and parameters recorded. The regime is varied
across datasets: 8 to 20 varieties, 40 to 90 concepts, root ages from 3000 to
8000 years, borrowing probabilities in {0, 0.15, 0.3}, and missingness fractions
in {0, 0.08, 0.2}. For each dataset two chains of the dated model are run and
summarized.

## Measures

Clade recovery is the fraction of true clades whose posterior support exceeds
0.5. Credible-region coverage is the fraction of datasets whose true root age
falls inside the 90 percent posterior interval; honest uncertainty means this is
near 0.9.

## Pilot results

A five-dataset pilot with short chains gives a mean clade recovery of 0.41 and a
90 percent root-age interval coverage of 0.8 (four of five). The moderate clade
recovery reflects the short pilot chains and the small, borrowing-perturbed trees;
coverage near the target is early evidence that the interval width is roughly
calibrated. The full study runs 200 datasets with the frozen budget; its summary
and per-dataset records are written to `runs/recovery/recovery.json`. Coverage
materially below 0.9 would indicate over-confident intervals and is reported as a
failure rather than hidden.

## What recovery cannot establish

Recovery under the model's own generative assumptions is a necessary check, not a
sufficient one: a posterior narrower than it should be can still be wrong under
misspecification. Chapter 12 therefore includes deliberate model violations, and
Chapter 9's posterior predictive checks test fit to features the generative model
does not directly target.
