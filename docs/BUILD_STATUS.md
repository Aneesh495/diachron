# Build status

Updated as the project progresses. Each gate resolves to `passed`, `failed`,
`blocked`, or `unverified` with a pointer to the evidence that justifies it.

## Active phase

Phase 1-2: repository foundation, Rust likelihood core, and tiny-state
correctness oracles.

## Gate table

| Gate | Description | Status | Evidence |
|---|---|---|---|
| G01 | Evidence profile (dimensions, licenses, attestations) | unverified | pending ingestion + synthetic corpus |
| G02 | Pruning matches enumeration on 1000+ tiny cases | passed | `tests/test_likelihood_oracle.py`: 1200 cases agree < 1e-9 |
| G03 | Sampler agrees with exact finite posteriors | passed | `tests/test_mcmc_exact.py`: discrete + Gamma + Jacobian control |
| G04 | Ascertainment correction validated | passed | `tests/test_ascertainment.py`: p(variable) matches sim, missing!=absent |
| G05 | 200+ known-truth recovery datasets | path verified (pilot) | `make recovery`: 5-dataset pilot clade-recovery 0.41, 90%CI coverage 0.8; full 200-run pending |
| G06 | Four-chain inference with convergence | unverified | `make sample` + `make diagnose` |
| G07 | Held-out concept/language prediction | path verified | `make predict` + ppc; p-values show tree fits sharing, model misses missingness |
| G08 | 12+ sensitivity configurations | path verified | `make sensitivity`: 12 configs; root age is calibration-dependent, stable to tree prior/ascertainment |
| G09 | Geographic honesty in UI/exports/thesis | passed (UI/export) | map view 3-kind legend (modern/attestation/inferred), exporter tags kinds; thesis pending |
| G10 | Borrowing recoverability limits | passed | `tests/test_geo_borrowing.py`: intensity stays low without signal, rises with real co-presence |
| G11 | Interactive workbench on real artifacts | passed | `ui/` renders 4 linked views on saved bundles; headless screenshot in docs/screenshots |
| G12 | Exact-state recovery, stale schema rejection | passed (resume) | `tests/test_recovery.py`: exact resume + corrupt-checkpoint reject |
| G13 | Semantic verifier negative controls | passed | `tests/test_verifier_controls.py`: 5 controls (fabricated likelihood, negative age, synthetic-as-observed, modern-as-ancient, corrupt digest) rejected by subprocess verifier for correct reason |
| G14 | Thesis claim-to-evidence closure | passed | `docs/thesis/`: 15 chapters + claims.json (9 claims, all resolve to evidence) |
| G15 | Clean-checkout release + hosted publication | unverified | release pipeline |

## Toolchain (verified this session)

- macOS arm64, 8 cores, 16 GB RAM
- Python 3.14.2 (venv), numpy 2.5, scipy 1.18, pandas 3.1, pyarrow 26, matplotlib 3.11
- Rust 1.88, cargo 1.88, maturin 1.16 (PyO3 bindings)
- Node 24.21, npm 11.19
- gh 2.101 authenticated as Aneesh495

## Last verified push

Repository: https://github.com/Aneesh495/diachron
Remote main tip verified at commit 30bcffd (Rust core + bindings).

## Progress log

- Phase 1-2 foundation: Rust core (expm, CTMC models, Felsenstein pruning,
  independent enumerator, serializable RNG, MCMC moves, geo kernels) compiles,
  21 unit tests pass, and builds as a CPython 3.14 extension. The pruning vs
  enumeration agreement is the start of G02; the Python oracle campaign will
  scale it to 1000+ cases.
