# Known-truth recovery

The current 200-dataset run used protocol 2 at executable commit
`8c47e96bab0f326d8476a90c418304018397df06`. All 200 datasets completed.
The central 90% root-age intervals contained the true age in 144 datasets
(72%). The mean fraction of nontrivial true clades with a sampling frequency
above one half was 0.4132. Root ages started from the derived calibration
midpoint, with an independently initialized topology in each chain.

The intervals undercovered in this workload. The run does not distinguish
short-chain mixing, simulator misspecification, and calibration misspecification
as explanations. The model treats the simulated data as a dated tree, even
when the simulator injected borrowing. This is a test of model robustness
under deliberate misspecification.

The earlier run used a different sampler and reported 94% coverage. Its
results predate the corrected topology proposal and observed-mask
ascertainment. Retain it as a historical result. Use protocol 2 for current
comparisons.

- [Run summary and stratification](recovery-v2.json)
- [Per-dataset ages and clade masses](recovery-v2-rows.jsonl.gz)
