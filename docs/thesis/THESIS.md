# Diachron: evidence and inference for language family history

## Abstract

Diachron implements pruning and posterior sampling for historical lexical
data. Its evidence model keeps lexical forms, source cognate judgments,
loan flags and uncertain dated observations distinct. The primary study
uses 152 IE-CoR varieties and 170 meanings. A bounded ABVD study covers
40 Austronesian varieties, and a separate EDH cohort contributes 175
individual Latin inscriptions.

The likelihood conditions on the actual character retention rule and
observed-tip masks. A finite contact mixture marginalizes one possible
copy event per character. The workbench links saved tree draws with
lexical observations and source records. Geographic root draws retain
their distinct label and show the diagnostic status of the underlying run.

The observed chains did not meet convergence thresholds. In the revised
200-dataset synthetic recovery study, nominal 90% root-age intervals
contained the true age in 72% of datasets and mean nontrivial clade
recovery was 41%. These results leave historical chronology and geographic
comparisons open.


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

Protocol 2 executed all 200 known-truth datasets at commit
`8c47e96bab0f326d8476a90c418304018397df06`. Each dataset was evaluated with
two chains, 400 warmup iterations and 1200 sampling iterations. Dates were
generated with three levels of attestation coverage. The sampler started from
a random tree and a calibration midpoint. The all-tip clade was excluded
from scoring.

The nominal 90% root-age intervals covered 144 of 200 true ages. Mean
nontrivial clade recovery was 0.4132, using a sampling frequency threshold
of 0.5. The interval result is below its nominal target. The observation model,
topology mixing and the simulator all contribute possible explanations. The
experiment does not isolate them.

The results directory includes raw sampled root ages and true-clade
frequencies per dataset, together with the seed, regime and source commit.
These permit direct recomputation of coverage and recovery without rerunning
inference. The earlier 94% result used a different sampler and selection
correction and remains available as an archived result.


# Chapter 9: Observed study results

## Evidence

The pinned IE-CoR primary study retains 152 varieties and 170 meanings. It has
24,989 cells with a lexical record, 851 missing lexical cells, and 2,255
loan-associated judgments. It contains 44 historical varieties with source
date distributions. The bounded ABVD study retains 40 Austronesian varieties
and 210 concepts. It has 7,355 lexical cells, including 1,169 with no cognate
judgment. These dimensions follow the frozen rules in the observed-data
protocol and do not depend on inferred trees.

The EDH cohort contributes 175 individual records with reported date intervals,
localities and bibliographies. The cohort covers Latin epigraphy and is indexed
to the IE-CoR Latin entry for inspection. It contributes no automatic
calibration. The IE-CoR dated models use the wordlist authors' terminal
variety distributions instead.

## Inference status

Each observed diagnostic used four independent starts, 400 warmup iterations,
1200 sampling iterations and thinning by ten. These checks test the whole path
through saved draws and do not set the convergence gate. Historical estimates
remain unverified. The source-defined historical tips are terminal leaves, so
the result also cannot adjudicate whether a recorded variety was an ancestor of
another sampled variety.

Held-out prediction and sensitivity operate as distinct experiments. Their
scores are not inferred from convergence or the retrospective fit of these
models.


# Chapter 10: A bounded contact model

Each lexical character begins as a binary pattern generated on the dated tree.
A categorical variable selects either the inherited pattern, with weight 1, or
one directed contact edge with weight equal to its nonnegative intensity. When
an edge is selected, the donor's latent tip state replaces the recipient's.
The candidate network is fixed by a modern-proximity screen. Directions are
hypotheses with separate intensities, independent exponential priors and equal
prior rates.

For an observed recipient, the copied pattern is scored by observing its
state at the donor and marginalizing the recipient's inherited state. If both
tips were observed and differ, the copy component has probability zero. If
the recipient is missing, copying changes no observed state. For the selection
event, an observed recipient adds one count to the donor and stops contributing
its own inherited count. The likelihood sums over contact and rate categories
before conditioning on retention. The no-contact model follows as all
intensities tend to zero.

The enumeration checks include missing donors and competition between two
edges. They validate the calculation, and the synthetic comparison tests
sensitivity to injected sharing on candidate edges. Independent historical
contact remains a separate question. The model permits only one terminal copy
per character, and a repeated match can arise from shared inheritance or a
copy. The model does not identify contact date, chains of transmission, or
contact outside the candidate set.

IE-CoR's loan classes and ABVD's form loan flags are auditable annotations.
The observed likelihood can be recomputed after removing those forms.
The contact component models visible cognate states and carries no likelihood
for the existence of an annotation.


# Chapter 11: Location evidence and inference

Diachron uses three distinct inputs or outputs. Contemporary lexical varieties
can carry catalog reference points. Independently sourced inscriptions carry
reported historical localities and dating uncertainty. The geographic model
produces latent samples conditional on its tree, rates, and modern reference
points. The workbench shows the categories separately. An ancient object may
have traveled before deposition, and its origin is not a boundary of a speech
community.

The planar contrasts model integrates internal locations by Gaussian message
passing and conditions on measured modern reference points. Missing tips
supply no geographic contrast. This approximation measures squared degrees
and is sensitive to the extent of the region. The region model instead uses a
symmetric distance-decayed transition generator and exact discrete pruning.
Unresolved region assignments can be marginalized.

The workbench's planar root cloud is drawn offline from the conditional root
distribution for each retained geographic tree and diffusion rate. It shows
points in a ±250-year band around the timeline position. Observations retain
their full date ranges. Root-only sampling is an explicit scope limit of the
current map export. A fitted regional model still needs its own offline root
state sampler.

None of the geographic diagnostics establish a homeland. Coordinate spread,
between-chain agreement, and stability to the reference-point model remain
empirical checks.


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
