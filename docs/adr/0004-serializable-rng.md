# ADR 0004: Carry a serializable RNG for exact resumption

## Status
Accepted.

## Context
Checkpoint and resume must reproduce a chain bit-for-bit, which requires capturing
and restoring the RNG state exactly. Rust's standard library has no stable,
serializable RNG, and NumPy's generator state is awkward to thread through the
language boundary on every draw.

## Decision
Implement a xoshiro256** generator with an explicit four-word state in Rust,
seeded via splitmix64, exposed to Python with `snapshot` and `from_state`. The
chain uses only this RNG for its stochastic decisions.

## Consequences
Exact resumption works: an interrupted-and-resumed chain is bit-identical to an
uninterrupted one (gate G12). The cost is that we maintain our own RNG rather than
reusing a library one, and the generator must be used consistently (no stray
NumPy calls in the chain's accept/reject path, or reproducibility would break).
