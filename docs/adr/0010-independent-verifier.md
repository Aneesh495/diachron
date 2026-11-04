# ADR 0010: A verifier separate from the result summarizer

## Status
Accepted.

## Context
Evidence bundles could be tampered with or could claim results their raw records
do not support. A verifier that trusted the summarizer's own pass booleans, or
shared its code, would not be an independent check.

## Decision
Implement the verifier (`diachron.verify`) separately from the summarizer. It
recomputes decisive quantities from raw records (draw digest, likelihood replay
through the saved engine, split-Rhat), checks geographic honesty, and rejects a
bundle for a semantic reason. It shares format schemas but not the decisive code
paths (it has its own Newick parser and Rhat implementation).

## Consequences
A tampered result is caught: altering a stored root age without updating the
manifest digest is rejected for a digest mismatch (a negative control in the
integration test). The cost is duplicated parsing and diagnostic code, which is
intentional; the duplication is the independence. The verifier never retrains,
regenerates, or trusts saved booleans.
