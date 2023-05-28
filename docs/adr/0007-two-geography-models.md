# ADR 0007: Offer two geography models and compare them

## Status
Accepted.

## Context
Geographic evidence varies in resolution. Forcing a single model risks either
pretending to resolve coordinates the data cannot support or discarding
continuous signal where it exists.

## Decision
Implement both a region CTMC (discrete latent regions, reusing the discrete
pruning path) and a Brownian-on-tree model (continuous coordinates, exact Gaussian
pruning). Compare their posterior predictions rather than declaring one correct.

## Consequences
The coarse-evidence case uses regions without overclaiming precision; the
continuous case uses diffusion where attestations constrain movement. The cost is
two models to maintain and a comparison to report rather than a single number. The
region model approximates its distance-decayed generator by an F81 chain to stay
on the fast pruning path, which is documented as an approximation and cross-checked
against the Brownian model rather than asserted as exact.
