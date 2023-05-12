# Diachron: posterior inference over language family history under uncertain evidence

## Abstract

Diachron is a research platform for studying how language families diversify,
borrow words, and shift their geographic footprint, built around original
probabilistic inference rather than narrative. The core computes a Felsenstein
pruning likelihood of cognate-class evolution under a continuous-time Markov
model, composes it with clock, calibration, tree-prior, geographic, and bounded
borrowing terms into an explicit posterior, and samples that posterior with an
original Metropolis-Hastings kernel whose moves carry derived and tested proposal
ratios and Jacobians.

The thesis makes two kinds of statement and keeps them apart. Mathematical
statements concern the algorithms: the pruning recursion computes the stated
marginal likelihood, and the sampler leaves the stated target invariant. These
are proved against independent tiny-state enumeration (agreement within 1e-9 over
1200 random cases) and against exact finite posteriors. Historical statements are
empirical claims with uncertainty: a posterior over a root age or a latent
homeland is a distribution with credible regions, prior sensitivity, and
counterevidence, never a recovered fact.

A central finding, established on labeled-synthetic known-truth corpora, is that
the inferred chronology is identifiable only through calibration: the root-age
estimate shifts materially when calibrations are removed or halved, while it is
stable to the tree-prior strength and the ascertainment correction. Posterior
predictive checks show the tree model reproduces the pairwise cognate-sharing
structure of the data but does not reproduce the missingness pattern, which the
likelihood handles but the generative check does not, a misspecification we report
rather than hide. Geographic conclusions are presented as posterior clouds; the
platform refuses to collapse a broad or multimodal location posterior to a single
homeland pin, and separates modern reference coordinates, dated attestations, and
inferred latent locations in every view and export.

The contribution is an auditable, reproducible inference pipeline with honest
uncertainty, not a claimed reconstruction of any particular family's history.


# Chapter 1: Research questions and hypotheses

These questions and hypotheses were frozen before the full experiments were run
(see `docs/EXPERIMENT_PLAN.md`). They separate what the evidence can identify
from what it cannot, and each resolves to a concrete comparison in a later
chapter.

## Predeclared questions

**Q1. Geographic versus purely phylogenetic predictive value.** Does conditioning
on geography improve held-out predictive density over a tree-only model, or does
the tree alone carry the signal? Identifiability concern: geographic tip
coordinates are modern reference points, so any geographic gain may reflect
present-day areal structure rather than history.

**Q2. Borrowing sensitivity.** When a known loan set is present, does a
borrowing-aware model recover the contact signal, and does removing loans from the
inherited characters change the tree? Limit: the borrowing model uses a bounded
candidate network, so it cannot recover unrestricted reticulation.

**Q3. Ascertainment bias.** Does the variable-character ascertainment correction
change inferences materially, and does switching missing-to-absent produce the
distinct likelihoods the model requires?

**Q4. Calibration dependence of chronology.** How much does the inferred root age
depend on the calibration priors versus the character likelihood? Hypothesis H4:
without calibration the time scale is unidentifiable and the root age is driven by
the prior.

**Q5. Rate heterogeneity.** Does a relaxed clock change divergence-time estimates
relative to a strict clock, and is the relaxed clock identifiable from this data?

**Q6. Prior sensitivity of geography.** Does the geographic prior (spherical
diffusion variance, or region granularity) drive the latent-location posterior
more than the data?

**Q7. Posterior multimodality.** Do independent chains agree on the tree
posterior, or is there multimodality that scalar diagnostics would miss?

**Q8. Influence of single observations.** Does removing one influential language
or one attestation move an apparently stable homeland or age estimate?

## Strongest initial hypothesis and its counterargument

The strongest initial hypothesis a reader might hold is that a coherent, tightly
dated homeland can be read off the posterior. Chapter 13 argues against this
directly: the platform's own sensitivity results show the chronology is
calibration-dependent and the location posterior is often broad, so a single
confident homeland claim is not supported by this evidence even when a point
estimate looks plausible. Whether that counterargument ultimately wins for a given
family is itself an empirical question the platform is built to ask, not assume.


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


# Chapter 3: Evidence, annotation, and the character transform

## Data contract

Lexical data enters as typed records: concepts, languages (with Glottolog codes
and modern reference coordinates), lexemes (with original script and
transliteration provenance), and cognate judgments (with confidence, loan and
partial-cognate flags, annotator, and source). Dated attestations carry a dating
interval in years before present, a source reference and page, an observation
type (attestation, inscription, reconstruction, trade), a location with
uncertainty, a language-assignment confidence, and explicit conflicting-record
links. Referential integrity and value ranges are checked on ingestion.

## CLDF ingestion

The CLDF reader locates the Wordlist tables by their component URL and resolves
each role-bearing column by its CLDF property, so it works across datasets that
name columns differently. It records the dataset license and DOI. Restricted
lexical collections are ingested from a directory acquired under their own
license and are never redistributed; the public release ships an acquisition
recipe and derived artifacts the license permits.

## The character transform and its audit trail

The transform from judgments to the character matrix is explicit and reversible
on paper. Each cognate class within a concept becomes one binary character whose
state for a language is present, absent, or missing by the rules in Chapter 4.
Every present and absent cell carries an evidence link back to the lexemes and
judgment sources that produced it, so an audit sample can be traced from a matrix
cell to the underlying annotation. Singleton classes are excluded and counted;
loans are recorded and, in the borrowing sensitivity analysis, removed from the
inherited characters to measure their effect.

## Evidence profile

The primary study corpus has 112 eligible varieties and 110 concepts with 11333
observed cells and 140 dated attestations across four observation types; the
second study has 36 varieties and 90 concepts with 60 attestations and its own
calibration regime. Both corpora are labeled synthetic known-truth: their history
is known, which is what the recovery and coverage studies need, and the synthetic
flag propagates into every artifact and view so they are never presented as
observed. The real-data path ingests CLDF collections of comparable scale when a
license permits.

## Correlated missingness and ascertainment

The transform preserves the information needed to test whether missingness,
replacement rate, and concept selection correlate with age or geographic
sampling. The corpora deliberately include uneven coverage and region-biased
attestation sampling so inference must reflect that bias rather than assume even
coverage.


# Chapter 4: Generative models

This chapter writes the complete generative model and its posterior, and shows
the limiting cases in which the full model reduces to a simpler one.

## Unit of inheritance

An observation attaches to a (language, concept) cell. The transform in
`diachron.io.characters` maps cognate judgments to binary presence/absence
characters: one character per cognate class within a concept. For language L and
class K the state is present (L has a lexeme judged in K), absent (L has lexical
data for the concept but no form in K), or missing (L has no lexical data for the
concept). Singletons (present in one language) are excluded as carrying no
internal-structure signal and counted in the ascertainment accounting. Loans are
annotated, not folded silently into inheritance; partial cognates are their own
classes.

## Character evolution

Each class evolves along the tree as a two-state continuous-time Markov chain with
instantaneous rate matrix

    Q = [ -gain   gain ]
        [  loss  -loss ]

so a lineage gains the class at rate `gain` and loses it at rate `loss`. The
stationary distribution is proportional to `(loss, gain)`. The transition matrix
over a branch of length t has the closed form

    P(t) = [ pi0 + pi1 e   pi1 - pi1 e ]
           [ pi0 - pi0 e   pi1 + pi0 e ],   e = exp(-(gain+loss) t),

with `pi0 = loss/(gain+loss)`, `pi1 = gain/(gain+loss)`. This is validated
against the matrix exponential of Q in the Rust test suite.

## Pruning likelihood

For one character the probability of the observed tip states, marginalized over
internal node states, is the Felsenstein pruning likelihood. Define the
conditional likelihood vector `L_v(s) = P(data below v | v in state s)`. Tips are
indicators of the observed state, or the all-ones vector when missing (which is
why missing and absent differ). At an internal node v with children c,

    L_v(s) = prod_c ( sum_j P_sj(t_c) L_c(j) ).

The character log likelihood is `log sum_s pi(s) L_root(s)`. Chapter 5 gives the
log-space scaling that keeps this stable on deep trees and the ascertainment
conditioning.

## Clock and calibration

A strict clock scales branch durations in years to expected substitutions by a
single rate; a relaxed uncorrelated-lognormal clock gives each branch its own
rate from a lognormal with shared mean and log-standard-deviation sigma. As sigma
goes to zero the relaxed clock reduces exactly to the strict clock. Calibrations
place densities on node ages: an aggregate root-age interval derived from the
oldest reliable attestation, which is the minimal structure that makes the strict
clock's time scale identifiable.

## Tree prior

Node times follow a constant-size coalescent: with k lineages the pairwise
coalescent rate is `k(k-1)/(2N)` and the waiting time to the next coalescence is
exponential with that rate. The density is written out in
`diachron.models.priors.coalescent_log_prior`.

## Geography

Lineage location is a latent state. Two models are compared rather than assumed. A
region CTMC moves a lineage among a small set of named regions. A Brownian-on-tree
model diffuses latitude and longitude with variance sigma^2 per year; its marginal
likelihood of the tip coordinates, integrating out internal locations, is the
Gaussian independent-contrasts pruning, the continuous analogue of the character
recursion. Tip coordinates enter as uncertain observations with their own error;
they are modern reference points, never historical boundaries.

## Borrowing

A bounded contact model places a borrowing intensity on each edge of a fixed
candidate network, with an exponential prior pulling intensities toward zero so
the model prefers no borrowing unless the data support it. This is a finite,
identifiable alternative to unrestricted reticulation, which the platform
deliberately does not claim.

## The full posterior

    log p(tree, times, rates, geo, borrow | data)
      = sum_characters log p(character | tree, times, rates)  [pruning]
        - n_variable * log p(variable | model)                 [ascertainment]
      + log p(node ages | calibrations)                        [calibration]
      + log p(tip coords | tree, times, sigma^2)               [geography]
      + log p(borrow intensities)                              [borrowing prior]
      + log p(tree)                                            [coalescent]
      + log p(rates)                                           [clock prior]
      + const.

## Limiting cases

- Borrowing intensities to zero: the borrowing model reduces to the dated model.
- Geographic sigma to infinity: locations become uninformative and the geography
  term is constant, recovering the dated model.
- Relaxed-clock sigma to zero: the relaxed clock reduces to the strict clock.
- No calibration: the calibration term drops and the time scale is identified only
  up to the clock-rate/time trade-off, which Chapter 12 shows empirically.


# Chapter 5: Likelihood derivations and numerical stability

## Pruning recursion correctness

Claim (mathematical). For a rooted tree with branch lengths, a reversible
two-state CTMC with transition matrices P(t), and observed tip states, the pruning
recursion computes exactly

    P(observed tips) = sum over all internal-state assignments of
                       pi(root state) * prod_edges P(parent state -> child state).

Proof sketch. The recursion is the distributive law applied to the sum over
internal states. Factor the joint probability as a product over edges; each
internal node's state appears only in the edges incident to it, so the sum over
that state can be pushed inward past all terms not involving it. Processing nodes
in post-order performs exactly these inner sums, yielding `L_v(s)` as the partial
sum over the subtree below v conditional on v's state. At the root the remaining
sum against pi is the full marginal. A missing tip contributes the all-ones vector
because its unobserved state is summed over all values, which is the marginal over
that tip; an observed tip contributes an indicator.

Verification (empirical). An independent enumerator
(`rust/diachron-core/src/enumerate.rs`, and a separate Python oracle in
`tests/oracle.py`) sums the joint probability over all internal-state assignments
with no shared code path. Over 1200 random trees with 2 to 5 tips, random rates,
and 15% missing tips, the pruning and enumeration log likelihoods agree to within
1e-9 (gate G02, `tests/test_likelihood_oracle.py`). The agreement of two
independently written routines in two languages is the evidence that the fast
algorithm is correct.

## Log-space scaling

Products of many small probabilities underflow. At each internal node the
conditional vector is divided by its maximum entry, and the log of that maximum is
accumulated into a per-node log scale; the final log likelihood adds the
accumulated scales to `log sum_s pi(s) L_root(s)`. This keeps the stored vector
O(1) regardless of tree depth. A caterpillar tree of 200 tips, which underflows
without scaling, returns a finite log likelihood with it
(`tests/` and `cli stress`).

## Missing versus absent

A missing tip marginalizes over its state; an absent tip pins the state to 0.
Because marginalizing sums the two indicator cases, the missing likelihood is
strictly larger than either pinned case. The Rust and Python tests both confirm
`L(missing) > L(absent)` for a fixed surrounding pattern (gate G04). This is the
mechanism by which "we have no data here" and "the class is absent here" produce
different likelihoods, as the evidence contract requires.

## Ascertainment conditioning

Only variable characters are kept. Conditioning each kept character on being
variable divides its likelihood by `p(variable) = 1 - p(constant)`, where
`p(constant)` is the summed probability of the all-absent and all-present patterns
under the same tree and model. We compute `p(constant)` by running the pruning
likelihood on the two constant patterns and summing. Validation: the analytic
`p(variable)` matches the empirical fraction of forward simulations that come out
variable, to within 0.005 over 400000 simulations (gate G04,
`tests/test_ascertainment.py`). The correction matches the selection process
actually applied; a correction for a different scheme would not reproduce the
simulated fraction.

## Geographic Gaussian pruning

The Brownian-on-tree likelihood uses independent contrasts: at each internal node
the children's location estimates are combined by precision weighting, the added
variance is the harmonic combination of the child branch variances, and each
child-pair contrast contributes a Gaussian density. This integrates out all
internal locations in closed form and gives the exact marginal likelihood of the
tip coordinates under Brownian motion, the continuous counterpart of the discrete
pruning recursion.


# Chapter 6: Sampler correctness

The sampler is a Metropolis-within-Gibbs cycle over scalar parameters and tree
structure. This chapter states each move's proposal ratio, derives the two that
carry a Jacobian, and reports the tests that validate them against exact targets.

## Acceptance rule

A proposed state is accepted with probability `min(1, exp(d))` where

    d = log p(proposed) - log p(current) + log_hastings,

and `log_hastings = log q(current|proposed) - log q(proposed|current) + log|J|`.
Non-finite proposals are rejected. The rule is validated on a finite five-state
target: sampled frequencies match the exact normalized probabilities to a maximum
absolute deviation below 0.01 over 400000 iterations (gate G03,
`tests/test_mcmc_exact.py`).

## Scalar scale move

For a positive scalar the move proposes `x' = x * exp(lambda (u - 1/2))`,
`u ~ U(0,1)`. The proposal density is symmetric in log space, so the Hastings
ratio is the Jacobian of the transform, `log(x'/x) = log(factor)`. This is tested
two ways. First, with the correct Jacobian the move targets a Gamma(3, 2) density
with sampled mean 6 and variance 12 (analytic values) within predeclared
tolerances. Second, a negative control runs the same sampler with the Jacobian
omitted and confirms the stationary mean is detectably wrong, so the passing test
is not vacuous.

## Node-time slide

An internal node's age is moved by a uniform step within a window, reflected to
stay strictly between its oldest child and its parent. Because the window is
recomputed consistently and out-of-window proposals are rejected, the move is
symmetric and the Hastings term is zero. Branch lengths are derived from node
ages, so time ordering is maintained automatically.

## Root scale move

All internal node ages are multiplied by a common factor `exp(lambda(u-1/2))`.
Scaling k independent coordinates contributes a Jacobian of `factor^k`, so the
Hastings term is `k log(factor)`. Invalid results (negative branch) are reverted
and rejected.

## Topology moves

**NNI.** A node's sibling is swapped with one of the node's children. Under
symmetric selection the move is its own inverse, so the Hastings term is zero. The
move is rejected if it would make a child older than its parent.

**Bounded SPR.** A node and its parent edge are pruned; the collapsed parent is
reinserted on a target edge whose time interval spans the pruned node's crown age.
Targets are restricted to edges outside the pruned subtree and excluding the
collapsed node and the sibling edge, which is what prevents cycles; the Hastings
ratio is the ratio of valid-target counts before and after. An earlier version
omitted these exclusions and could create a cycle, which made likelihood traversal
blow up; the fix and a 3000-move integrity test that checks for cycles, a single
root, a stable node count, and valid times are in `tests/test_moves.py`.

## Checkpointing and exact resumption

The chain's RNG is a serializable xoshiro256** with a four-word state. A
checkpoint writes the RNG words, tree, parameters, frozen step sizes, iteration
counters, and retained draws atomically via a temp-file rename. Resuming
reconstructs an identical chain. A chain interrupted and resumed partway produces
a trajectory bit-identical to an uninterrupted run: same final log posterior to
1e-9, same topology, same draw count, and no duplicated retained draw (gate G12,
`tests/test_recovery.py`). Adaptation of step sizes runs only during warmup and is
frozen before any retained sample, so the retained chain targets the stationary
posterior.


# Chapter 7: Experimental protocol

The protocol is frozen in `docs/EXPERIMENT_PLAN.md`. This chapter summarizes it.

## Models compared

Four observation models score the same character matrix with the same
ascertainment correction, so differences in fit come from added structure, not
from changing the data: tree-only, dated, geography-conditioned, and
borrowing-aware.

## Chains and stopping

The full inference campaign runs at least four independently initialized chains
per model. Each chain starts from its own random tree and its own RNG seed.
Warmup adapts step sizes; adaptation is frozen before any retained draw. Stopping
budgets are frozen after a pilot that measures per-step cost and acceptance.

## Convergence criteria

Declared scalar endpoints are root age, gain, loss, and clock rate. A model is
reported converged only when split-Rhat is below 1.01 and bulk and tail ESS
exceed 400 for every endpoint, alongside topology mixing summaries and
between-chain clade agreement. Short diagnostic runs that fail these thresholds
are reported as not converged; their draws are never presented as a converged
posterior.

## Evaluation units

Known-truth recovery uses independent synthetic datasets as the evaluation unit,
spanning tree shape, missingness, borrowing, and dating regimes. Held-out
prediction holds out whole concepts so no cognate class spans the train/test
partition. Masks are applied before any fitting.

## Diagnostic versus acceptance profiles

A CPU diagnostic profile uses short chains for quick feedback; the full
acceptance profile uses long chains and writes a sealed evidence directory. The
two are named distinctly and a diagnostic run cannot set the acceptance flag.

## Decision rules

A conclusion is reported as stable only if it holds across the predeclared
sensitivity configurations; otherwise it is reported as prior- or
calibration-dependent. Failed and inconclusive analyses are retained in the index
of negative results, not dropped.


# Chapter 8: Synthetic recovery

Known-truth recovery answers two questions: can the sampler recover a tree it
generated under its own model, and are its credible regions honest?

## Design

Each synthetic dataset is drawn with a coalescent tree, cognate-class evolution
under the inference model, optional borrowing, Brownian geography, and dated
attestations, with the true tree and parameters recorded. The regime is varied
across datasets: 8 to 20 varieties, 40 to 90 concepts, root ages from 3000 to
8000 years, borrowing probabilities in {0, 0.15, 0.3}, and missingness fractions
in {0, 0.08, 0.2}. For each dataset two chains of the dated model are run and
summarized.

## Measures

Clade recovery is the fraction of true clades whose posterior support exceeds
0.5. Credible-region coverage is the fraction of datasets whose true root age
falls inside the 90 percent posterior interval; honest uncertainty means this is
near 0.9.

## Pilot results

A five-dataset pilot with short chains gives a mean clade recovery of 0.41 and a
90 percent root-age interval coverage of 0.8 (four of five). The moderate clade
recovery reflects the short pilot chains and the small, borrowing-perturbed trees;
coverage near the target is early evidence that the interval width is roughly
calibrated. The full study runs 200 datasets with the frozen budget; its summary
and per-dataset records are written to `runs/recovery/recovery.json`. Coverage
materially below 0.9 would indicate over-confident intervals and is reported as a
failure rather than hidden.

## What recovery cannot establish

Recovery under the model's own generative assumptions is a necessary check, not a
sufficient one: a posterior narrower than it should be can still be wrong under
misspecification. Chapter 12 therefore includes deliberate model violations, and
Chapter 9's posterior predictive checks test fit to features the generative model
does not directly target.


# Chapter 9: Results on the study corpus

Results here are on the primary study corpus (112 varieties, 110 concepts, 1576
variable characters, 140 attestations). The corpus is labeled synthetic
known-truth; the same pipeline runs on CLDF collections ingested under their own
licenses.

## Inference and convergence

A diagnostic four-chain run of the dated model (1200 retained iterations per
chain after 400 warmup) produces a root-age posterior with mean about 7240 years
and active topology exploration: 60 to 72 distinct topologies per chain, a
topology switch rate near 0.55, and a between-chain clade agreement of 0.76. The
scalar endpoints on this short run have split-Rhat between 1.1 and 2.8 and ESS in
the tens, so the run is reported as not converged; it is a diagnostic profile, not
the acceptance campaign. The honest reading is that a short run of a 112-tip tree
space has not mixed, and more iterations are required before the posterior is
reported as converged.

## Held-out prediction

Held-out-concept prediction holds out whole concepts so no cognate class spans
the partition, fits on the training characters only, and scores the held-out
characters by posterior predictive log density per character. The comparison
between the dated and tree-only models, and the leakage guard, are recorded in
`runs/predict/heldout.json`. The design, not a single number, answers Q1:
geographic or dated structure earns its place only if it improves held-out
density over the tree-only baseline.

## Posterior predictive checks

Simulating characters from posterior draws and comparing summary statistics to
the observed data gives a Bayesian p-value near 0.18 for mean pairwise cognate
sharing, so the tree model reproduces who-shares-cognates well. The missingness
rate has a p-value of 0 because the generative check does not model the
missingness process that the likelihood handles by marginalization; this is a
reported misspecification of the generative check, not of the likelihood. The
mean class size is slightly over-predicted (p near 0.93). These results are the
kind of honest mixed fit the platform is built to surface.


# Chapter 10: Borrowing and contact

## Model and its limits

The borrowing model places a borrowing intensity on each edge of a fixed
candidate contact network, with an exponential prior pulling intensities toward
zero. This is a bounded, identifiable model. It deliberately does not attempt
unrestricted network reconstruction: with a fixed candidate set it can ask whether
a specific posited contact improves the fit, but it cannot discover arbitrary
reticulation, and the platform does not market it as doing so.

## What it can and cannot recover

On known-truth data where a loan was injected between geographic neighbors, the
borrowing-aware model can detect elevated sharing on candidate edges that the tree
alone does not explain, provided the true contact lies on or near a candidate
edge. When the true contact is absent from the candidate set, the model cannot
recover it; the correct outcome in that case is that the borrowing channel stays
near its zero-intensity prior and the conflict shows up instead as reduced tree
support, not as a false contact claim.

## Loan removal sensitivity

Removing loan-annotated forms from the inherited characters and refitting measures
how much the tree depends on borrowed material. This is one of the predeclared
sensitivity configurations (Chapter 12). A tree that changes substantially under
loan removal is flagged as borrowing-sensitive, which is a finding about
identifiability, not a defect.

## No unsupported causal history

The platform reports borrowing as a posterior over intensities on named candidate
edges with stated priors, never as a narrative of which people borrowed from whom.
A detected intensity is evidence that a candidate contact improves fit under the
model; it is not a proof of historical contact, and Chapter 13 keeps that
distinction explicit.


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


# Chapter 12: Sensitivity and robustness

Twelve predeclared configurations each change one lever and refit, so stable
findings can be separated from prior- and calibration-dependent ones. The levers
are: coalescent population size at three values (tree-prior strength),
ascertainment on versus off, singleton inclusion, loan removal, a half-calibration
subset, no calibration, and the four observation models.

## Headline result: chronology is calibration-dependent

The root-age estimate shifts materially under the calibration-related levers
(halving the calibrations, removing calibration, and the tree-only model) and
under the geography model, while it is stable to the tree-prior strength and to
the ascertainment correction. This confirms hypothesis H4: without calibration the
time scale is not identifiable and the chronology is driven by the prior. The
practical consequence is that any dated conclusion from this platform must be
reported together with its calibration, and a chronology presented without stating
its calibration dependence would be misleading.

## Ascertainment and tree-prior stability

That the root age is stable to turning ascertainment off and to a factor-of-nine
change in the coalescent population size is itself a finding: the character
likelihood, not those priors, is driving the age estimate given a fixed
calibration. Stability here does not mean the ascertainment correction is
unnecessary, only that it does not dominate the dated age; it does change
likelihood values and the relative fit of constant-heavy versus variable-heavy
regimes, which matters for model comparison.

## Deliberate misspecification

Beyond the twelve levers, recovery includes datasets generated with borrowing and
missingness regimes the fitted model does not match, so the study tests whether a
confidently narrow posterior survives model violation. A posterior that narrows
below its honest width under misspecification is reported as over-confident.

## Reporting rule

A conclusion appears in the summary as stable only if it holds across the relevant
levers. Calibration-dependent and geography-dependent conclusions are labeled as
such. The full configuration table is written to
`runs/sensitivity/sensitivity.json`.


# Chapter 13: Counterarguments

This chapter argues against the strongest initial hypothesis a reader might bring:
that Diachron recovers a confident, tightly dated homeland for a family.

## The chronology argument against confident dating

The sensitivity results (Chapter 12) show the inferred root age is
calibration-dependent: it moves materially when calibrations are weakened or
removed. A confident date therefore reflects the calibration as much as the
lexical data. Unless the calibration itself is defended from primary evidence, a
narrow dated claim overstates what the cognate data alone identify. This argument
does not say dating is impossible; it says a date is only as good as its stated
calibration, and the platform forces that statement into the open.

## The homeland argument against a single pin

Geographic tip coordinates are modern reference points. A latent-location
posterior inferred from them can be coherent and still be driven by present-day
areal structure rather than by history. When the posterior is broad or multimodal,
reporting a single homeland pin discards exactly the uncertainty that matters. The
influence check compounds this: if dropping one language or attestation moves the
estimate, the pin was never stable. The platform's refusal to collapse the cloud
is the methodological form of this counterargument.

## The mixing argument against trusting short runs

A short run of a large tree space can look settled on scalar summaries while the
topology has not mixed. The diagnostics (Chapter 9) show exactly this on the
diagnostic-profile run: scalar Rhat and ESS fail the thresholds and the run is
reported as not converged. A result read off an unconverged chain is not a result,
and the acceptance gate refuses to certify one.

## When the counterargument might lose

For a family with strong, independently defended calibrations, dense attestation,
and chains that pass the convergence thresholds across independent initializations,
a dated posterior with a genuinely narrow credible interval and a concentrated,
removal-stable location posterior would warrant a confident claim. Whether that
holds is an empirical question per family; the platform is built to test it, not
to assume the answer. The honest default, absent that evidence, is the
counterargument.


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
