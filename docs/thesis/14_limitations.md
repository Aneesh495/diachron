# Chapter 14: Limitations and what the evidence cannot identify

## Modeling limitations

The borrowing model uses a fixed candidate network; it cannot discover arbitrary
reticulation and does not try to. The geographic observation model treats tip
coordinates as noisy present-day references; it cannot turn modern points into
ancient boundaries, and a geographic gain in predictive density may reflect areal
rather than historical structure. The character model collapses a concept's forms
into binary presence/absence classes under a stated transform; polymorphism,
partial cognacy, and synonymy are represented but not fully modeled as a joint
process.

## Inference limitations

Topology mixing on large trees is the binding constraint. The diagnostic runs
here did not converge, and convergence on 100-plus-tip trees requires long chains
and careful proposal tuning. More samples do not rescue a stuck topology chain;
the platform reports non-convergence rather than presenting unconverged draws.

## Identifiability limits

The chronology is not identifiable without calibration. The clock-rate and
time-scale trade off, and the data fix only relative times until a calibration
anchors them. A broad or multimodal location posterior is a genuine limit of the
evidence, not a rendering artifact. Marginal-likelihood model comparison is
offered only where the estimator is validated; raw likelihoods from different
observation models are never compared as if they were Bayes factors.

## Data limitations

The headline corpora are labeled synthetic known-truth, which is ideal for
recovery and coverage but is not a substitute for observed data; the real-data
path depends on CLDF collections the user acquires under their own licenses.
Attestation records are sparse and biased by design, and the inference reflects
that sparsity rather than smoothing it away.

## What a reader should not conclude

A reader should not conclude from this work that any particular family's tree,
date, or homeland is established. The established results are mathematical: the
likelihood is correct against independent enumeration, the sampler targets exact
finite posteriors, the ascertainment correction matches its selection process,
and interrupted chains resume exactly. The historical results are posteriors with
stated uncertainty, counterevidence, and the explicit calibration and mixing
caveats above.
