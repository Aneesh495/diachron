# ADR 0011: The workbench reads saved artifacts and never infers

## Status
Accepted.

## Context
An interactive workbench that recomputed inference on interaction (for example,
re-sampling when the timeline moves) would be slow, non-reproducible, and could
trigger unbounded work from a user gesture.

## Decision
The workbench loads a single exported JSON bundle and reads already-sampled draws.
Moving the timeline, selecting a clade, or switching models only reads this data.
The Python `export-workbench` command writes the bundle from saved posterior
artifacts.

## Consequences
Interactions are cheap and bounded, and what the user sees is exactly the saved
posterior, so the UI cannot drift from the evidence. The cost is an export step
and a potentially large JSON file (a few megabytes for two models on the primary
corpus), which is acceptable for a research tool and is capped (the matrix is
limited to 400 columns for responsive rendering). The honesty labels (synthetic
banner, geographic legend) live at this boundary so every view inherits them.
