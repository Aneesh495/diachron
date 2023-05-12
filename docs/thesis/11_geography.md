# Chapter 11: Geography

## Observation model

A catalog coordinate is a present-day reference location, not a historical
homeland or boundary. An ancient attestation may be an isolated inscription,
elite writing, or a trade document rather than a whole-population distribution.
Diachron encodes the observation type and treats tip coordinates as uncertain
observations with their own error, so the geographic likelihood never assumes a
modern point marks an ancient boundary.

## Two models compared

The region CTMC moves a lineage among coarse named regions; it is appropriate
when the evidence cannot resolve coordinates and avoids over-claiming precision.
The Brownian-on-tree model diffuses coordinates continuously and integrates out
internal locations exactly by Gaussian pruning. Rather than declare one map
historically correct because it looks coherent, the platform compares their
posterior predictions; the region model is preferred where the data are coarse,
the diffusion model where attestations constrain continuous movement.

## Posterior clouds, not pins

A latent location is reported as a posterior cloud or credible region. When the
posterior is broad or multimodal, it stays a cloud; the platform refuses to
collapse it to a single homeland pin. The workbench draws the inferred location
as a density of sampled points whose spread grows with inferred temporal depth,
and the map legend keeps modern reference points, dated attestations, and
inferred latent locations visually distinct at every time step (gate G09).

## Sensitivity and influence

The geographic prior (diffusion variance, region granularity) is one of the
predeclared sensitivity levers; a location posterior driven more by the prior than
the data is reported as such. The influence of a single language or attestation is
tested by removal: an apparently stable homeland estimate that moves when one
influential observation is dropped is not a stable finding. This directly serves
Q6 and Q8.
