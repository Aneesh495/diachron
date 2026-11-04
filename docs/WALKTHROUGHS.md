# Source-linked walkthroughs

Four complete execution paths traced through the source, including an error and
recovery path. File references point at the modules that implement each step.

## Walkthrough 1: a character to a likelihood (the happy path)

1. `python/diachron/io/simulate.py::simulate_dataset` produces a `Dataset` of
   typed records (or `io/cldf.py::load_cldf_wordlist` for real CLDF data).
2. `io/characters.py::build_character_matrix` transforms cognate judgments into
   binary characters, resolving present/absent/missing and attaching an
   `EvidenceLink` per cell.
3. `infer/campaign.py::prepare_study` aligns character columns to the tree's
   taxon order and derives calibrations, tip coordinates, and borrow candidates.
4. `models/posterior.py::PosteriorModel.character_log_likelihood` calls the Rust
   `block_log_likelihood` (via `rust/diachron-core/src/bindings.rs`) over all
   columns, then adds the ascertainment correction from
   `models/ascertainment.py`.
5. The Rust `pruning.rs::log_likelihood_block` runs the Felsenstein recursion
   with log-space scaling and returns the summed log likelihood.

## Walkthrough 2: one MCMC step

1. `infer/chain.py::Chain._step` picks a move via the project RNG
   (`rust/diachron-core/src/rng.rs`).
2. For a topology move, `infer/moves.py::MoveSet.spr` prunes and regrafts the
   tree, returning the Hastings term as the log ratio of valid-target counts.
3. `Chain._eval` calls `PosteriorModel.log_posterior`, which sums the character
   likelihood, calibration, geography, borrowing, tree prior, and clock prior.
4. `Chain._accept` applies the Metropolis rule; on rejection the tree arrays are
   restored from the pre-move snapshot.
5. `Chain._maybe_sample_hyperparameters` occasionally samples geographic sigma
   and borrowing intensities when those models are active.

## Walkthrough 3: a campaign to a workbench

1. `cli/sampling.py::run_sample` builds the primary corpus, prepares the study,
   and calls `infer/campaign.py::run_campaign` for four independently initialized
   chains.
2. `infer/campaign.py::campaign_bundle` packages the draws into a
   `PosteriorBundle` with a draw digest and the source commit.
3. `infer/artifacts.py::PosteriorBundle.save` writes `bundle.json` and
   `MANIFEST.json`.
4. `cli/workbench.py::export_workbench` reads the saved bundles and writes
   `workbench.json` with modern points, attestations, and the matrix, tagging
   each geographic kind.
5. `ui/src/main.ts` loads the JSON and mounts the four linked views; selecting a
   clade in `ui/src/views/treeView.ts` propagates through `ui/src/store.ts` to
   the map, timeline, and matrix.

## Walkthrough 4: an error and recovery path

1. A chain is interrupted mid-run. Its last atomic checkpoint
   (`infer/chain.py::Chain.checkpoint`) holds the RNG words, tree, parameters,
   frozen step sizes, and retained draws.
2. `Chain.resume` reconstructs an identical chain. If the checkpoint is missing
   its RNG state or has a stale schema, `resume` raises rather than silently
   proceeding (tested in `tests/test_recovery.py::test_corrupted_checkpoint_rejected`).
3. On a successful resume, continuing the chain reproduces the exact trajectory
   it would have taken uninterrupted
   (`tests/test_recovery.py::test_exact_resume`).
4. When the completed evidence is later verified, `verify/verifier.py` recomputes
   the draw digest and replays likelihoods; a tampered or inconsistent bundle is
   rejected for a semantic reason, with the five negative controls in
   `tests/test_verifier_controls.py` demonstrating the rejection path.
