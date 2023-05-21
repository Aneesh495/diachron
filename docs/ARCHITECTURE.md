# Architecture

Diachron is a layered system: a Rust numerical core, a Python orchestration and
modeling layer, a TypeScript workbench, and a documentation and evidence layer.
Each layer has a narrow contract with the ones it touches.

## Layers

### Rust core (`rust/diachron-core`)

The compute kernel, exposed to Python via PyO3. It owns the numerics that must be
fast and exact:

- `linalg`: dense small-matrix operations and a scaling-and-squaring Pade matrix
  exponential, validated against closed forms and scipy.
- `subst`: continuous-time Markov substitution models (binary gain/loss, F81
  k-state) with closed-form transition matrices.
- `tree`: the arena tree representation and post-order traversal.
- `pruning`: the Felsenstein pruning likelihood with log-space scaling and the
  missing-versus-absent distinction.
- `enumerate`: an independent brute-force likelihood used only to validate
  `pruning`.
- `rng`: a serializable xoshiro256** generator for reproducible, resumable chains.
- `mcmc`: proposal moves and the acceptance rule (also mirrored in Python for the
  driver).
- `geo`: haversine distance and the Brownian-on-sphere log density.
- `bindings`: the PyO3 surface; flat arrays cross the boundary, no opaque handles.

### Python layer (`python/diachron`)

- `io`: typed records, CLDF ingestion, the auditable character transform, the
  known-truth simulator, and the reference corpora.
- `models`: the timed tree, clocks, explicit prior densities, the ascertainment
  correction, and the composed `PosteriorModel` across four configurations.
- `infer`: the MCMC moves, the chain driver with checkpointing and exact resume,
  posterior artifacts, and the multi-chain campaign.
- `diagnostics`: split-Rhat, rank-normalized Rhat, Geyer ESS, bulk/tail ESS,
  MCSE, and topology mixing.
- `geo`: the region CTMC model and contact classification.
- `ppc`, `sensitivity`: posterior predictive checks, held-out prediction, the
  recovery study, and the predeclared sensitivity configurations.
- `registry`: the alternative-history registry and attestation inspector.
- `verify`: the independent semantic verifier.
- `cli`: the research command line that wires everything into runnable commands.

### TypeScript workbench (`ui`)

Reads saved posterior bundles and renders four linked views (tree, map, timeline,
lexical matrix) with selection propagation and uncertainty legends. It performs
no inference; the Python `export-workbench` command writes the bundle it reads.

### Documentation and evidence (`docs`, evidence directories)

The specification, this architecture, ADRs, research notes, diagrams, the thesis,
and the experiment plan. Evidence directories hold sealed campaign outputs that
the verifier checks independently.

## Key data path: a character to a likelihood

1. `io.cldf` or `io.simulate` produces a `Dataset` of typed records.
2. `io.characters.build_character_matrix` transforms cognate judgments into binary
   characters with present/absent/missing states and an evidence trail.
3. `infer.campaign.prepare_study` aligns character columns to the tree's taxon
   order and derives calibrations and tip coordinates.
4. `models.PosteriorModel.character_log_likelihood` calls the Rust
   `block_log_likelihood` over all columns and adds the ascertainment correction.
5. `infer.chain.Chain` proposes moves, evaluates the composed posterior, accepts
   or rejects, and records thinned draws.
6. `infer.artifacts.PosteriorBundle` saves the draws with a digest and source
   commit; `cli.workbench.export_workbench` turns bundles into the UI's JSON.

## Ownership and concurrency

The Rust core is stateless across calls except for the `PyRng` object, whose state
the Python chain owns and checkpoints. Chains are run sequentially in the campaign
driver; there is no shared mutable state between chains, so independent
initialization is genuine. The UI is single-threaded and reactive through the
store.

## Failure and recovery

Chains checkpoint atomically (temp file plus rename) with the full RNG and tree
state. Resume reconstructs an identical chain. Corrupt or schema-mismatched
checkpoints are rejected rather than silently accepted. The verifier recomputes
decisive quantities from raw records and rejects tampered or inconsistent evidence
for a semantic reason.
