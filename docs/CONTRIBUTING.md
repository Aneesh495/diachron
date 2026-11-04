# Contributing

Diachron is a research platform; contributions should preserve its two
invariants: the numerical core stays correct and auditable, and the distinction
between observation, inference, and illustration stays visible.

## Setup

```bash
make bootstrap
make build
make test
```

## Expectations for a change

- New likelihood or sampler behavior needs a correctness test. For the engine,
  prefer a differential check against an independent implementation (the
  enumeration oracle, an exact finite target, or a closed form) over a test that
  re-derives the same algorithm.
- A topology or structure-changing move needs an integrity test (no cycles, one
  root, stable node count, valid times) in addition to its Hastings-ratio
  derivation.
- Any historical claim added to the thesis needs a claim-to-evidence entry in
  `docs/thesis/claims.json`.
- Public prose follows the writing rules: explain mechanisms and measured
  behavior, no em dashes, no hype, no boasting about size or counts.

## Running the gates

- `make likelihood-check` runs the pruning-vs-enumeration oracle (G02).
- `make test` runs the Rust unit tests and the full Python suite, including the
  exact-MCMC, ascertainment, resume, and verifier-control tests.
- `make sample` then `make diagnose` runs a campaign and reports convergence;
  `make acceptance` and `make verify` run and check the full campaign.

## Code layout

See `docs/ARCHITECTURE.md`. The Rust core is in `rust/diachron-core`, the Python
layer in `python/diachron`, and the workbench in `ui`.

## What not to do

- Do not add padding, repeated adapters, or generated variation to inflate size.
- Do not present synthetic data as observed, collapse a broad location posterior
  to a single pin, or compare raw likelihoods from different observation models
  as if they were Bayes factors.
- Do not weaken a convergence threshold or shorten a required workload to force a
  gate to pass; report non-convergence instead.
