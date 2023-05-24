# Note 9: Geography

Geography is the part of the platform most at risk of overclaiming, so its design
is mostly about what it refuses to do.

It refuses to treat a modern catalog coordinate as an ancient boundary. Tip
coordinates enter as uncertain present-day references, and the observation type of
an attestation (inscription, reconstruction, trade, or general attestation) is
recorded because an isolated inscription is not a population centroid. It refuses
to apply a single historical map; the region model and the Brownian model are
compared rather than one being declared correct because it looks coherent. And it
refuses to collapse a broad or multimodal location posterior to a single homeland
pin: the workbench draws the inferred location as a density whose spread grows
with inferred temporal depth, and the legend keeps modern, attested, and inferred
marks distinct.

The region model treats location as a discrete latent state over coarse named
regions and reuses the discrete pruning likelihood, which is the right granularity
when the evidence cannot resolve coordinates. The Brownian model diffuses
coordinates continuously and integrates out internal locations exactly by Gaussian
pruning. A reported result is the comparison of their posterior predictions, not a
single map.

The influence check is the sharpest honesty tool: an apparently stable homeland
that moves when one influential language or attestation is removed was never
stable. A geographic conclusion survives this note only if it survives that
removal, and the platform is built to run it rather than to assume the estimate is
robust.
