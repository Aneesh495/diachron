# Execute and audit the full profile

Start from a clean checkout at the published source commit. Install the pinned
Python, Rust and UI tools using the README. Use a host with space for the saved
posterior bundles and checkpoints. The full profile runs all declared fits and
records each command exit. Unconverged fits retain their results and fail the
corresponding scientific checks.

```sh
.venv/bin/diachron acquire --out runs/inputs --study iecor --study abvd
make acceptance RUN_ID=full-observed PROFILE=full ACCEPTANCE_ARGS='--inputs runs/inputs --jobs 4'
make verify EVIDENCE=evidence/full-observed VERIFY_ARGS='--expected-commit SOURCE_COMMIT'
```

Pass `--release-policy` with the operator's separately maintained policy when
checking release readiness. The acceptance record keeps the policy outcome.
The exact policy bytes remain in the operator's storage.

The full profile records the frozen input reports, executable source inventory,
platform, command exits, per-model draws, final checkpoints, browser trace,
wheel and per-gate decisions. Failed commands retain stdout and stderr. Run
`release --evidence evidence/full-observed` after independent verification.

To inspect existing executions, pass `--evidence-plan` in place of `--inputs`.
An audit uses recorded evidence and leaves full execution required for release.
The plan is a JSON object with `schema` equal to 1 and `source_commit` equal to
the clean verification checkout. Paths resolve relative to the plan file.
Use these sections for the available evidence.

| Section | Fields |
| --- | --- |
| `inputs` | `iecor`, `abvd`, `catalog`, `edh_source`, `edh_geography` |
| `inference`, `prediction`, `sensitivity` | `directory`, `source_commit` |
| Archived `sensitivity` | `archive_directory`, `asset_manifest`, `source_commit` |
| `recovery` | `summary`, `rows`, `source_commit` |
| `workbench` | `html`, `qa` |

The verifier checks the pinned bytes, recomputes recovery summaries and
convergence, replays held-out probes and runs the semantic controls. It reports
missing components and mismatched sources separately. A report of incomplete
research returns a nonzero status and preserves the gate details.

For archived sensitivity, use a trusted asset manifest with `source_commit`
and an `assets` object mapping archive filename to SHA-256. The verifier
restores and checks one setting at a time. The archive inventory and the
posterior diagnostics remain separate checks.
