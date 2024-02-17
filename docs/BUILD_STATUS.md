# Build status

The code, observed data pipeline, and diagnostic workflows are runnable.
Full acceptance remains blocked. All four full-budget observed fits completed and missed the convergence criteria.

| Gate | Status | Evidence and open work |
| --- | --- | --- |
| G01 Evidence | passed | The audit rebuilt 152 primary varieties, 170 meanings, 40 secondary varieties, 206 source-proposed cognacy links and 175 individual EDH records. The records cover Latin inscriptions and carry assignment limitations. See [Observed data](OBSERVED_DATA.md). |
| G02 Likelihood | passed | The tiny-tree oracle tests cover 1,200 independently enumerated cases. Rust tests also check large-tree scaling and a wide polytomy. |
| G03 Sampler | passed | Independent enumeration checks all four composed finite targets. Production acceptance frequencies, reversible topology moves and scale Jacobians meet the declared tests. See [Sampler correctness](thesis/06_sampler_correctness.md). |
| G04 Ascertainment | passed | Exact selection on observed masks with and without singleton filtering is checked against enumeration. Contact copies and rate mixtures are included in the selection normalization. |
| G05 Recovery | failed | The earlier and matched studies each completed 200 datasets. Nominal 90% root intervals covered 72% and 88.5% of true ages, respectively. See [Recovery](results/recovery_200.md). |
| G06 Inference | failed | All 16 full-budget chains completed and replayed. Each model missed scalar and clade agreement thresholds. See [full inference](results/full-inference.md). |
| G07 Prediction | unverified | The lexical split precedes character construction, and the language split fixes its evaluation vocabulary from training data. The observed campaign is running. |
| G08 Robustness | unverified | Twenty observed configurations include alternate cognacy, catalog constraints, calendar priors, clocks and geography. The observed campaign is running. |
| G09 Geographic provenance | passed | Map marks distinguish modern points, sourced historical localities and model-sampled root locations. Unknown locations are kept off the map. |
| G10 Borrowing | passed | The bounded contact model matches independent enumeration with missing donors and two candidate edges. The [Brittonic fruit audit](notes/11_historical_case_studies.md) exposes a missed candidate and cross-class transfer stays outside this model. |
| G11 Workbench | passed | The four-model workbench passed [offline desktop and mobile checks](results/workbench-full-browser.json), including the five-step evidence tour and catalog comparisons. |
| G12 Resume | passed | Checkpoints reject mismatched source and model identities. Exact resume and equality of serial and parallel draws have regression tests. |
| G13 Verifier | passed | The verifier accepted the original full dated bundle and [rejected four rehashed semantic edits](results/observed-semantic-controls.json). The audit independently reproduced the gate decisions and retained the research failures. |
| G14 Thesis | passed | Fifteen Markdown chapters and a standalone HTML renderer passed the claim-ledger audit. Empirical results distinguish verified execution from convergence and calibration. |
| G15 Release | blocked | The public source and bounded demonstration can be verified from a clean checkout. The release command remains gated on independently verified full acceptance. |

## Continue

- Use [Observed data](OBSERVED_DATA.md) to reconstruct source inputs and
  diagnostic runs from the pinned files.
- Run `make test`, `make likelihood-check` and `make demo` for the bounded
  regression suite.
- Preserve the existing observed runs before changing initialization or
  convergence protocol. Compare a revised sampler against the original
  protocol and record failures at chain level.
