# Diagrams

Editable Mermaid sources (`.mmd`) with checked-in rendered SVG exports. Re-render
with `make docs` (requires the mermaid CLI; see `docs.py`). Each diagram is kept
synchronized with the source it describes.

| File | Covers |
|---|---|
| 01_system_context | How a researcher, data sources, and Diachron interact |
| 02_components | Component map across the Rust core, Python layer, and UI |
| 03_evidence_lineage | Raw lexeme judgment to model contribution (one data path) |
| 04_character_evolution | The two-state CTMC for a cognate class (domain math) |
| 05_pruning_messages | Conditional-likelihood message flow up the tree |
| 06_proposal_moves | MCMC moves and their Hastings terms |
| 07_calendar_constraints | Dates, calibration, and clock relationships |
| 08_borrowing_network | Bounded candidate contact network with contact kinds |
| 09_posterior_map_render | Three kinds of geographic mark kept distinct |
| 10_chain_lifecycle | Warmup, adaptation freeze, sampling, checkpoint (state lifecycle) |
| 11_failure_recovery | Checkpoint, resume, rejection, and verification paths |
| 12_testing_and_acceptance | How gates map to checks and evidence |
