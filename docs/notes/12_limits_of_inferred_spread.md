# Note 12: Limits of inferred spread

This note states plainly what the geographic inference cannot establish, so that a
coherent-looking spread animation is never mistaken for an observation.

An inferred latent location is a posterior over where a lineage might have been,
conditioned on modern reference coordinates and a diffusion or region model. It is
not a record of where speakers were. Three gaps separate it from history. First,
the tip coordinates are modern: a diffusion fit to them can recover present-day
areal structure and project it backward, which looks like history but is not.
Second, the attestations that could anchor the past are sparse and biased toward
documented, often elite or written, contexts; they do not sample the speaker
distribution. Third, the model has no access to barriers, migrations, or contacts
it was not given, so a smooth inferred spread reflects the smoothness of the
diffusion prior as much as any evidence.

The platform therefore reports spread as a posterior cloud with its uncertainty,
labels illustrative animations as illustrative, and keeps them visually separate
from attested observations at every time step. It runs the single-observation
influence check and reports when the spread estimate is not robust to dropping one
language or attestation.

The honest bottom line is that Diachron can say which spread hypotheses are
consistent with the data and how uncertain each is; it cannot, from cognate data
and modern coordinates, certify a particular migration. Presenting one as certain
would be exactly the error the platform is built to prevent.
