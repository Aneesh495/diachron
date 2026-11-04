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
| G05 | 200+ known-truth recovery datasets | unverified | `make recovery` |
| G06 | Four-chain inference with convergence | unverified | `make sample` + `make diagnose` |
| G07 | Held-out concept/language prediction | unverified | `make predict` |
| G08 | 12+ sensitivity configurations | unverified | `make sensitivity` |
| G09 | Geographic honesty in UI/exports/thesis | unverified | UI + export labels |
| G10 | Borrowing recoverability limits | unverified | `tests/test_borrowing.py` |
| G11 | Interactive workbench on real artifacts | unverified | `ui/` + serve path |
| G12 | Exact-state recovery, stale schema rejection | passed (resume) | `tests/test_recovery.py`: exact resume + corrupt-checkpoint reject |
| G13 | Semantic verifier negative controls | unverified | `diachron/verify` |
| G14 | Thesis claim-to-evidence closure | unverified | `docs/thesis/` |
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
