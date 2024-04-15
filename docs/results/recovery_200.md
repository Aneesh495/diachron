# Known-truth recovery

The 200-dataset misspecification run used protocol 2 at executable commit
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
ascertainment. Retain it as a historical result.

- [Run summary and stratification](recovery-v2.json)
- [Per-dataset ages and clade masses](recovery-v2-rows.jsonl.gz)

## Matched observation model

Protocol 3 ran at `2df7fe5a3607fe671af075e3084a1a809dc6e652`. It generated
serial trees, dated tips and retained columns from the fitted model. Each
of the 200 datasets used four independent chains. Nominal 90% intervals
covered 177 true root ages (88.5%). Mean nontrivial clade recovery was 24.7%.
The independently recomputed root rank histogram had 61 of 200 cases in the
two highest bins. The raw ranks use serially dependent draws.

| Regime | Datasets | Root coverage | Mean clade recovery |
| --- | ---: | ---: | ---: |
| Vertical | 101 | 86.1% | 25.1% |
| Candidate copy | 99 | 90.9% | 24.2% |

The overall result missed the declared 90% target. The simulator, calibration,
tree starts and chain budgets all changed between protocols. The difference
in coverage cannot isolate any one change.

- [Run summary](recovery-v3-matched.json)
- [Per-dataset truth and sampled root ages](recovery-v3-matched-rows.jsonl.gz)
- [Independent coverage, rank and generator replay](recovery-v3-matched-verification.json)

```sh
.venv/bin/diachron verify-recovery \
  --summary docs/results/recovery-v3-matched.json \
  --rows docs/results/recovery-v3-matched-rows.jsonl.gz \
  --expected-commit 2df7fe5a3607fe671af075e3084a1a809dc6e652
```

## Protocol 4 paired matched campaign

The same 200 simulated truths, masks and seed rules were retained. Each fit
used four chains, 1,600 warmup iterations and 6,400 sampling iterations, with
every tenth draw saved. The joint time and clock block was enabled.

| Measurement | Protocol 3 | Protocol 4 |
| --- | ---: | ---: |
| Root coverage at nominal 90% | 177/200 (88.5%) | 170/200 (85%) |
| Mean clade recovery | 24.7% | 28.8% |
| Median root bulk ESS | 20.7 | 176.7 |
| Median root rank R-hat | 1.149 | 1.021 |
| Median root interval width in years | 9,262 | 9,760 |

The exact saved truths agree across all pairs. Coverage held in 168 cases,
failed in both fits in 21, fell out of the interval in nine, and moved into
it in two. Interval calibration remains below the target. The sampler and
budget changed together, and the paired comparison cannot isolate either
effect.

- [Run summary](recovery-v4-matched.json)
- [Per-case truth, root draws and clade membership](recovery-v4-matched-rows.jsonl.gz)
- [Independent replay](recovery-v4-matched-verification.json)

```sh
.venv/bin/diachron verify-recovery \
  --summary docs/results/recovery-v4-matched.json \
  --rows docs/results/recovery-v4-matched-rows.jsonl.gz \
  --expected-commit a135e990262e3d2224b1cd9ec0be8a6dac2831a3
```
