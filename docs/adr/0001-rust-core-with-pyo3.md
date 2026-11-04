# ADR 0001: Rust numerical core behind a Python orchestration layer

## Status
Accepted.

## Context
The inner loops (pruning likelihood, MCMC moves, matrix exponential) run millions
of times and must be fast and exact. The orchestration (data ingestion, model
composition, diagnostics, CLI) benefits from Python's ecosystem and iteration
speed.

## Decision
Implement the numerical kernel in Rust and expose it to Python via PyO3, with a
flat array interface. Keep the engine original rather than wrapping an existing
phylogenetics library.

## Consequences
The hot path is fast (about ten thousand character likelihoods per second on a
60-tip tree) and auditable, since every step is our code. The cost is a language
boundary crossed on every MCMC step and a build dependency on a Rust toolchain
and a PyO3 version that supports the target Python. PyO3 0.22 did not support
Python 3.14, which forced an upgrade to 0.25 and a feature-gated extension-module
flag so plain `cargo test` can still link without a Python interpreter.
