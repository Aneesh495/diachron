# Diachron 0.1.0 preview

Download `Diachron_Preview.html` and open it in a browser. The workbench runs
offline and displays saved draws from the bounded synthetic demo. The corpus
has 16 varieties, 80 concepts, and 40 dated attestations. Use the timeline and
attestation inspector to see how dates and source records affect the display.
[Preview downloads](https://github.com/Aneesh495/diachron/releases/tag/v0.1.0-preview)

`Diachron_Thesis.html` describes the models and links to the frozen experiments.
`diachron-preview-evidence.tar.gz` contains the demo summary, replay bundle,
and verification records. `SHA256SUMS` identifies the downloadable assets.

## Replay the saved run

From a checkout of the preview tag, build the extension and extract the
downloaded evidence archive to `runs/preview-evidence`. The replay checks the
bundle against the checked out source identity and recomputes saved scores.

```sh
git checkout v0.1.0-preview
make bootstrap
make build
mkdir -p runs/preview-evidence
tar -xzf /path/to/diachron-preview-evidence.tar.gz -C runs/preview-evidence
make verify EVIDENCE=runs/preview-evidence/demo VERIFY_ARGS="--expected-commit $(git rev-parse HEAD)"
```

## Research status

The preview demonstrates the pipeline from synthetic data through sampling,
score replay, and offline exploration. Observed inference, prediction, and
sensitivity missed their frozen convergence thresholds. The latest matched
recovery covered 85% of true root ages against a nominal 90% target. See the
gate table for evidence and limits on interpreting the historical results.
[Build status](BUILD_STATUS.md)
