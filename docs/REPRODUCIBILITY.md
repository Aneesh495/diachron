# Reproducibility

Everything in Diachron is reproducible from a clean checkout with a pinned
toolchain. This document lists the exact commands and the seeds and profiles that
make runs repeatable.

## Toolchain

- Python 3.11 or newer (developed on 3.14).
- Rust stable (developed on 1.88).
- Node 18 or newer (developed on 24) for the workbench.

## From a clean checkout

```bash
make bootstrap     # create .venv and install pinned Python tooling
make build         # compile the Rust core (maturin develop --release)
make test          # Rust unit tests + full Python suite
make demo          # bounded end-to-end demonstration (~8s)
```

## Determinism

- MCMC chains are seeded. The campaign uses base seed 1000 with a fixed
  per-chain multiplier; recovery datasets use 20000 plus the dataset index. All
  seeds are recorded in the posterior bundle and the evidence directory.
- The RNG is a serializable xoshiro256**; a chain resumed from a checkpoint
  reproduces the exact trajectory it would have taken uninterrupted.
- The synthetic corpora are produced by `numpy.random.default_rng` with fixed
  seeds (primary 20240601, second 20240602), so the data regenerates identically.

## Profiles

- CPU diagnostic: short chains for quick feedback; cannot set the acceptance flag.
- Full acceptance: `make acceptance RUN_ID=<id> PROFILE=full` runs the declared
  campaign into `evidence/<id>` and records `SOURCE_INPUTS.json`,
  `DATASETS.json`, `ENVIRONMENT.json`, `COMMANDS.jsonl`, `ACCEPTANCE.json`, and
  `MANIFEST.json`.

## Independent verification

```bash
make verify EVIDENCE=evidence/<id>
```

The verifier recomputes decisive quantities (draw digests, likelihood replay,
diagnostics) from the raw records and rejects a bundle that does not hold up. It
is implemented separately from the summarizer that wrote the evidence.

## Artifact formats

- Posterior bundle: `bundle.json` (model, dataset identity, per-chain draws) plus
  `MANIFEST.json` (draw digest and source commit).
- Workbench data: `workbench.json`, written by `diachron export-workbench`.
- Evidence directory: the JSON records listed under Profiles, plus per-model
  bundles and checkpoints.

## Recovery runbook

If a chain is interrupted, its last checkpoint is in the run's `checkpoints`
directory. Resuming is automatic in the campaign driver; a corrupt or
schema-mismatched checkpoint is rejected rather than silently accepted. See
`docs/notes/07_convergence.md` for how non-convergence is reported.
