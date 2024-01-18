# Full-budget observed inference

Protocol 3 ran four models on pinned IE-CoR v1.2.1. Each model retained four
independently seeded chains with 3000 warmup steps, 12000 sampling steps and
thinning by ten. All 16 chains completed. The executable source was
`01459af9ab79bca1c97ac2490624223a294292a2`.

| Model | Root split R-hat | Root bulk ESS | Root tail ESS | Maximum clade difference | Sampling seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| tree_only | 1.5265 | 7.6 | 24.0 | 1.00 | 675.0 |
| dated | 1.9876 | 5.9 | 21.8 | 1.00 | 318.2 |
| geography | 1.6234 | 6.3 | 10.7 | 1.00 | 412.9 |
| borrowing | 1.7859 | 5.9 | 14.2 | 1.00 | 1565.3 |

All four fits failed the frozen convergence thresholds. Each had at least
one clade with a frequency difference of 1.00 between chains. Root-age bulk
ESS stayed below eight. These results leave the observed inference gate open.
The diagnostic draws are available for examining chain behavior.

The independent verifier rebuilt the observed lexical columns and source
attestations. It checked the source inventory, complete budgets, final
checkpoints and 16 posterior scores per model. Those checks passed. The four
failed checks were the four convergence conditions. The
[machine-readable report](observed-full-01459af.json) records the checks,
source identities, elapsed sampling time and per-model draw digests.
