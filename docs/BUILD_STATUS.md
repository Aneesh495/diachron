# Build status

The code, observed data pipeline, and diagnostic workflows are runnable.
Full acceptance remains blocked. The reported sampling frequencies are from
diagnostic budgets, and the observed chains missed the convergence criteria.

| Gate | Status | Evidence and open work |
| --- | --- | --- |
| G01 Evidence | unverified | The observed dimensions and 175 individual records are documented in [Observed data](OBSERVED_DATA.md). Broader independent attestation coverage and adjudication of conflicting alternatives remain open. |
| G02 Likelihood | passed | The tiny-tree oracle tests cover 1,200 independently enumerated cases. Rust tests also check large-tree scaling and a wide polytomy. |
| G03 Sampler | unverified | Finite-target, scale Jacobian, and reversible topology move tests pass. Complete validation of the composed dated, geographic, and contact samplers remains open. |
| G04 Ascertainment | passed | Exact selection on observed masks with and without singleton filtering is checked against enumeration. Contact copies and rate mixtures are included in the selection normalization. |
| G05 Recovery | failed | All 200 revised datasets ran. Nominal 90% root intervals covered 72% of true ages. Mean nontrivial clade recovery was 41%. See [Recovery](results/recovery_200.md). |
| G06 Inference | failed | The observed diagnostic chains missed the split R-hat and ESS thresholds. Full converged model comparisons remain open. See [diagnostic data](results/observed-diagnostic.json). |
| G07 Prediction | unverified | The lexical split precedes character construction, and the language split fixes its evaluation vocabulary from training data. The observed comparisons and decision-grade estimates remain open. |
| G08 Robustness | unverified | Fifteen distinct sensitivity levers include calibration subsets, strict and relaxed clocks, site rate variation, two region resolutions and contact priors. Full observed sensitivity remains open. |
| G09 Geographic provenance | passed | Map marks distinguish modern points, sourced historical localities and model-sampled root locations. Unknown locations are kept off the map. |
| G10 Borrowing | unverified | Contact scoring matches independent enumeration with missing donors and two candidate edges. The [Brittonic fruit audit](notes/11_historical_case_studies.md) exposes a missed contact candidate. Cross-class comparative recovery remains open. |
| G11 Workbench | unverified | The saved-draw workbench passed [desktop and mobile checks](results/workbench-browser.json). The full research tour, taxonomic comparisons and general ancestral location views remain open. |
| G12 Resume | passed | Checkpoints reject mismatched source and model identities. Exact resume and equality of serial and parallel draws have regression tests. |
| G13 Verifier | unverified | Four rehashed semantic edits are rejected. The observed verifier can rebuild lexical columns, dates and inscriptions from pinned inputs. A complete campaign verifier remains open. |
| G14 Thesis | unverified | Fifteen Markdown chapters and a standalone HTML renderer are available. The empirical argument depends on the unfinished model comparisons. |
| G15 Release | blocked | The public source and bounded demonstration can be verified from a clean checkout. The release command remains gated on independently verified full acceptance. |

## Continue

- Use [Observed data](OBSERVED_DATA.md) to reconstruct source inputs and
  diagnostic runs from the pinned files.
- Run `make test`, `make likelihood-check` and `make demo` for the bounded
  regression suite.
- Preserve the existing observed runs before changing initialization or
  convergence protocol. Compare a revised sampler against the original
  protocol and record failures at chain level.
