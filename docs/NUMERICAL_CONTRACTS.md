# Numerical contracts

The precisions, units, and tolerances the engine promises, in one place.

## Units

- Node ages (`height`): years before present, nonnegative, tips usually at 0.
- Branch length: years (parent age minus child age).
- Substitution rates (`gain`, `loss`): expected events per year.
- Clock rate: substitutions per year; expected substitutions on a branch equal
  rate times branch length in years.
- Geographic diffusion `sigma2`: degrees squared per year (Brownian) or an
  implicit per-year rate (region CTMC).
- Distances: kilometers (great-circle via haversine, mean Earth radius).

## Precision

- All likelihoods and densities are double precision natural logs.
- The matrix exponential uses scaling-and-squaring with a degree-6 Pade
  approximant; on the small, well-conditioned generators Diachron builds it is
  accurate to near machine precision, cross-checked against closed forms and
  scipy.
- The pruning likelihood uses log-space scaling (per-node max factored out into
  an accumulated log scale) so it stays finite on deep trees.

## Tolerances (predeclared)

- Pruning vs enumeration (G02): max absolute log-likelihood difference below
  1e-9 over at least 1000 tiny cases.
- Sampler finite-target frequency (G03): max absolute deviation below 0.01.
- Scale-move Gamma recovery (G03): mean within 0.15, variance within 1.2 of
  analytic.
- Ascertainment (G04): analytic p(variable) within 0.005 of the simulated
  fraction over 400000 draws.
- Exact resume (G12): resumed log posterior matches uninterrupted within 1e-9;
  topology and draw count identical.

## Invariants

- A rate matrix has rows summing to zero; its exponential has nonnegative rows
  summing to one.
- A missing tip contributes the all-ones conditional vector; an absent tip
  contributes an indicator at state 0.
- A valid time tree has every child no older than its parent (zero-length
  branches allowed, negatives rejected).
- The RNG state is four 64-bit words; identical state yields identical streams.

## Error handling

- Non-finite proposed log posteriors are rejected, not propagated.
- Invalid trees (negative branch, cycle) are rejected by the move (which reverts)
  and would be caught by the integrity tests.
- Corrupt or schema-mismatched checkpoints raise on resume rather than loading
  partial state.
