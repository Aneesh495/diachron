# Chapter 2: Related work

Diachron builds on a body of phylogenetic and computational-linguistics method,
and it reimplements the algorithmic core rather than wrapping an existing engine
so that every step of the inference is auditable. The ideas below are attributed;
no published result is reproduced or claimed as Diachron's own.

## Phylogenetic likelihood

The pruning recursion for computing the likelihood of discrete characters on a
tree is due to Felsenstein. Diachron implements the two-state and k-state cases
with log-space scaling and an explicit missing-versus-absent distinction, and
validates the implementation against independent enumeration rather than against
another phylogenetics package.

## Bayesian phylogenetics and dating

Bayesian inference of dated trees with molecular clocks, node calibrations, and
relaxed-clock rate variation is standard in the phylogenetics literature and in
packages such as the BEAST family and MrBayes. Diachron does not use those
engines; it implements its own Metropolis-Hastings moves with derived proposal
ratios so the sampler is testable against exact finite targets.

## Convergence diagnostics

The convergence diagnostics follow the rank-normalized, split-chain R-hat and
effective-sample-size definitions of Vehtari, Gelman, Simpson, Carpenter, and
Burkner (2021). Diachron implements these directly and cross-checks them against
hand-computed small cases and analytic autoregressive examples.

## Computational historical linguistics

Cognate-based phylogenetic studies of Indo-European, including the dated analyses
associated with Gray, Atkinson, Bouckaert, and more recently Heggarty and
colleagues, motivate the data model (cognate classes per concept, dated
attestations as calibrations) and the questions about chronology and homeland.
Diachron treats those studies as sources of competing hypotheses and data
structure, not as ground truth; its real-data path ingests CLDF collections such
as IE-CoR under their own licenses rather than redistributing them.

## Data standards

CLDF provides the lexical data contract; Glottolog provides versioned language
identifiers and reference coordinates. Diachron pins dataset versions and DOIs and
keeps per-dataset licenses, because a shared schema does not grant one license to
every collection.

## What Diachron adds

The contribution is not a new algorithm but an integrated, honest, reproducible
pipeline: original likelihood and sampler with correctness evidence, an explicit
and auditable character transform, a bounded borrowing model with stated limits, a
geography layer that keeps observation types distinct, and a workbench and thesis
that never present inference as observation.
