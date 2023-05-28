# ADR 0002: Validate the likelihood against an independent enumerator

## Status
Accepted.

## Context
The pruning recursion is subtle and everything downstream depends on it. A test
that re-derives the same algorithm would not catch a shared conceptual error.

## Decision
Compute the character likelihood a second way, by brute-force enumeration over all
internal-node states, in both Rust (`enumerate.rs`) and a separate Python oracle
(`tests/oracle.py`), sharing no code with the recursion. Require agreement to
1e-9 over at least 1000 random tiny trees (gate G02).

## Consequences
A real differential test: a bug would have to be independently reproduced in an
enumeration to escape detection, which is far less likely than a transcribed
algorithm sharing a bug. The cost is that enumeration is exponential in internal
nodes, so the oracle is limited to tiny trees; this is acceptable because
correctness bugs show up on tiny trees, and separate stress tests cover large-tree
numerical stability.
