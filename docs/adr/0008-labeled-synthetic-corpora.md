# ADR 0008: Ship labeled-synthetic corpora, ingest real data by recipe

## Status
Accepted.

## Context
The restricted real lexical collections (IE-CoR, Lexibank) cannot be
redistributed. The platform still needs data at the required scale, and the
recovery and coverage gates need known ground truth.

## Decision
Ship two reproducible labeled-synthetic corpora at the required scale, with the
synthetic flag propagated into every artifact and view. Provide a schema-driven
CLDF reader that ingests real collections from a directory the user acquires under
the collection's own license, plus an acquisition recipe in the manifest.

## Consequences
The evidence profile is met with data that can be redistributed, and known truth
powers the recovery study. Synthetic data is never presented as observed history,
and the verifier rejects a bundle that claims a synthetic corpus is observed. The
limitation is that the headline numbers are on synthetic data; the real-data path
exists and works but its results depend on the user's licensed acquisition, which
the public repository cannot include.
