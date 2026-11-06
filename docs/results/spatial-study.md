# Predicting held-out locations on known trees

The executable source was `108a924cae9f3c6ad28fbcde345e24afdcfae105`.
The experiment fixed the tree in both models and predicted the same six
location categories. Each truth regime had 40 training trees for choosing
the rate multiplier and 120 test trees with three hidden tips each.

| Truth | Region CTMC negative log score | Spherical mesh negative log score | Paired difference and 95% interval |
| --- | ---: | ---: | ---: |
| Spherical diffusion approximation | 0.952 | 0.848 | −0.103 [−0.150, −0.057] |
| Six-region jump process | 0.382 | 0.415 | 0.033 [−0.001, 0.067] |

The mesh improved held-out prediction under spherical diffusion. The paired
interval under regional truth spans zero. Under spherical truth, nominal 90%
root-region sets covered 94.2% with the mesh and 99.2% with the regional
model. Their mean set sizes were 3.65 and 4.74 of six categories. The known
trees remove topology error from this experiment. Observed historical
predictions remain conditional on uncertain trees and their location inputs.

- [Aggregate scores and training scales](spatial-study.json)
- [All 240 held-out cases](spatial-study-cases.jsonl)
- [Independent score and matrix-exponential replay](spatial-study-verification.json)

```sh
.venv/bin/python scripts/verify_spatial_study.py \
  --summary docs/results/spatial-study.json --cases docs/results/spatial-study-cases.jsonl \
  --expected-commit 108a924cae9f3c6ad28fbcde345e24afdcfae105
```
