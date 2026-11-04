# Recovery study: 200 known-truth datasets (gate G05)

Run: `make recovery` with 200 datasets, two chains each, dated model. All datasets
are labeled synthetic known-truth. Raw per-dataset records are in
`runs/recovery_full/recovery.json`.

## Headline

- Root-age 90 percent credible-interval coverage: 0.94 (188 of 200). Target 0.90.
- Mean clade recovery (fraction of true clades with posterior support above 0.5):
  0.333; median 0.286.

## Coverage by borrowing regime

| borrow probability | datasets | root age in 90% CI | coverage |
|---|---|---|---|
| 0.0 | 67 | 64 | 0.96 |
| 0.15 | 70 | 65 | 0.93 |
| 0.3 | 63 | 59 | 0.94 |

## Reading

The coverage is at or slightly above the 0.90 target and holds across all three
borrowing regimes, including the ones that misspecify the fitted dated model. That
is the result the study is for: the root-age credible intervals are honestly
calibrated, and they do not become over-confident when the data were generated
with borrowing the model does not include.

Clade recovery is modest. Two reasons, both honest rather than defects to hide.
First, the per-dataset chains are short (1200 retained iterations), and topology
on even small trees needs more mixing, consistent with the convergence note.
Second, the borrowing-perturbed datasets inject shared cognates that conflict with
the tree, which correctly lowers support for the affected clades rather than
forcing a confident wrong tree. A recovery number inflated by longer chains would
be reported separately, with its own budget; this study reports what the frozen
budget produced.

## What this establishes and what it does not

It establishes that, under the model's generative assumptions and across a span
of borrowing and missingness regimes, the uncertainty Diachron reports for the
root age is calibrated. It does not establish clade recovery at production
quality, which needs longer chains, nor does it substitute for fit to real data,
which the posterior predictive checks probe separately.
