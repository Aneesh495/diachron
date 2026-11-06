# Held-out prediction and sensitivity

The observed campaign used pinned IE-CoR inputs at executable source
`d272a49aa4f57ee0d062b9fef0954c7e3e71759b`. Each of the 28 fits had four
chains, 3,000 warmup steps, 12,000 sampling steps and thinning by ten.
The prediction job took 7,128 seconds for both splits. The archive contains
the frozen partitions, scores, trajectories and checkpoints.
[Hosted execution](https://github.com/Aneesh495/diachron/actions/runs/42146898983)
[Raw archives](https://github.com/Aneesh495/diachron/releases/tag/research-42146898983-1)
[Archive digests](observed-research-assets.json)

## Prediction

The concept split withheld 34 concepts and scored 510 characters. The language
split withheld 30 varieties and scored 2,209 characters from the training-defined
vocabulary. Changes below are equal-weight means over concepts relative to the
tree-only fit. Positive values indicate higher scores on the saved draws.
Intervals resample concepts and reflect variation in their score differences.

| Model | Held-out concept change | Held-out language change |
| --- | ---: | ---: |
| Dated | −0.104 [−0.184, −0.015] | 0.102 [0.047, 0.150] |
| Geography | −0.004 [−0.071, 0.061] | 0.297 [0.236, 0.354] |
| Borrowing | −0.137 [−0.214, −0.059] | −0.032 [−0.082, 0.016] |

The local audit reconstructed both evaluation vocabularies and replayed
scores and final checkpoints for all eight fits. Every fit failed convergence.
The geographic scores differ by evaluation design. The values above summarize
diagnostic trajectories and the concept bootstrap omits uncertainty from
poor posterior mixing. [Full score report](observed-prediction-d272.json)

## Sensitivity

All 20 settings completed and passed archive, source, checkpoint and posterior
replay checks. Every setting failed convergence. Eight shifted the saved
root-age mean by more than 15% of the baseline mean. These included removing
ascertainment, halving or broadening date information, removing calibration,
planar and three-region geography, a tighter contact prior and a wider
diffusion prior. The undated setting's time scale is prior-dependent.
These differences describe the computed trajectories. Historical robustness
remains unresolved. [Sensitivity replay and per-setting results](observed-sensitivity-d272.json)

## Verify the archives

Use the pinned inputs from the [acquisition guide](../OBSERVED_DATA.md).

```sh
gh release download research-42146898983-1 --repo Aneesh495/diachron \
  --pattern '*.tar.gz' --dir runs/research-archives
.venv/bin/diachron restore-evidence \
  --archive runs/research-archives/prediction.tar.gz --out runs/restored-prediction \
  --sha256 6b7f28a8a99df9ba3d828d4f80be9d412c53fafe12d1759dd254a3de8dd637a7
.venv/bin/diachron verify --evidence runs/restored-prediction/fit \
  --expected-commit d272a49aa4f57ee0d062b9fef0954c7e3e71759b \
  --source runs/inputs/iecor --catalog runs/inputs/glottolog53.zip \
  --edh-source runs/inputs/edh --edh-geography runs/inputs/edh-geography.json
.venv/bin/diachron verify-sensitivity-archives --archive-dir runs/research-archives \
  --asset-manifest docs/results/observed-research-assets.json \
  --expected-commit d272a49aa4f57ee0d062b9fef0954c7e3e71759b \
  --source runs/inputs/iecor --catalog runs/inputs/glottolog53.zip \
  --edh-source runs/inputs/edh --edh-geography runs/inputs/edh-geography.json
```

The sensitivity command restores one archive at a time. Its nonzero exit
retains the scientific convergence failure even when artifact checks pass.
