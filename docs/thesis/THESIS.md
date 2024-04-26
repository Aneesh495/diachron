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

The observed chains missed convergence thresholds at full budget. Nominal
90% root-age intervals covered 72% of the 200 misspecified synthetic cases
and 88.5% of the 200 matched cases. Both held-out experiments and all 20
observed sensitivity settings completed. Their fitted draws leave historical
chronology, geographic prediction and robustness unresolved.


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

## Character and rate model

The unit is the membership of a language variety in a cognate class for a
particular concept. One source judgment can make the character present.
Absence requires full classification of the recorded forms in the cell.
A missing lexical cell or a cell with inadequate classification contributes
an unknown state. Loan judgments and multiple or partial classifications
remain in the evidence tables. They do not define mutually exclusive
states for a single meaning.

Let g and l be the gain and loss rates, and c the positive clock
multiplier. The two-state generator has off-diagonals g and l and
stationary distribution (l, g)/(g+l). Branch v has elapsed time
d_v = h_parent - h_v. The transition is exp(Q c d_v r_v a), where r_v
is one for the strict clock and a is one for the homogeneous site model.
Under the relaxed clock, log(r_v) follows a Normal with mean -s²/2 and
standard deviation s for fixed model dispersion s. Under the site
mixture, a is integrated over equal-probability discrete-gamma categories
scaled to arithmetic mean one. The clock multiplier and g and l enter
the likelihood as products. Proper priors regularize that redundancy.

## Observation and posterior

For character j, let x_j be its observed pattern and m_j its
observed-tip mask. Let A_j be the event of at least one observed absence
and either one or two observed presences, depending on the singleton
policy. For a fixed tree, define P_j as the mixture over rate categories
and possible contact edges. A contact component copies one donor's
latent terminal state to the recipient. It is averaged with the inherited
component before dividing by selection probability.

    log target = sum_j [log P_j(x_j | m_j) - log P_j(A_j | m_j)]
               + log p(ordered labeled tree | observed tip ages, N)
               + sum of dated tip and node calibration log densities
               + geographic log contrast or region likelihood
               + log p(g, l, c, r, geographic rate, contact intensities).

The denominator is omitted in the declared uncorrected sensitivity run.
Lognormal priors on g and l have median 0.001 and log standard deviation
1.5. The clock multiplier has median one and log standard deviation 0.5.
The planar geographic diffusion variance has median 0.004 squared degrees
per year and log standard deviation two. Contact weights have independent
exponential priors, with default rate ten. Each prior applies in the model
where its parameter is sampled.

For each waiting interval with k extant backward lineages, the constant-size
coalescent contributes survival under total hazard k(k-1)/(2N). Each labeled
merge contributes its own hazard 1/N. Historical tip dates introduce lineage
sampling events. The IE-CoR model uses terminal-tip densities in years before
2000 CE. The synthetic model can use an aggregate root-age calibration.
The ABVD tree-only analysis does not assign independent historical dates.

The dated baseline omits geography and contact. The geographic variant adds
planar contrasts or a symmetric region CTMC. The borrowing variant includes
contact copying. The contact likelihood approaches the dated lexical
likelihood as intensities approach zero. The relaxed clock approaches the
strict clock as its dispersion approaches zero.


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

Products of many small probabilities underflow. After each child contribution the
conditional vector is divided by its maximum entry, and the log of that maximum is
accumulated into a per-node log scale; the final log likelihood adds the
accumulated scales to `log sum_s pi(s) L_root(s)`. This keeps the stored vector
O(1) regardless of tree depth. A caterpillar tree of 200 tips, which underflows
without scaling, returns a finite log likelihood with it
(`tests/` and `cli stress`).

## Missing versus absent

A missing tip marginalizes over its state; an absent tip pins the state to 0.
Because marginalizing sums the two indicator cases, the missing likelihood is
at least as large as either pinned case. Equality is possible if the other state has zero mass. The Rust and Python tests both confirm
`L(missing) > L(absent)` for a fixed surrounding pattern (gate G04). This is the
mechanism by which "we have no data here" and "the class is absent here" produce
different likelihoods, as the evidence contract requires.

## Ascertainment conditioning

Condition on retention given the observed mask. The default rule needs
at least two observed presences and at least one observed absence.
When singletons are allowed, the positive threshold is one. Unobserved
tips contribute zero to both counts.

For each ancestral state, the normalization pass tracks six probabilities.
A bucket records the presence count saturated at two and a flag for any
absence. The transition step propagates child buckets through the CTMC.
The merge convolves saturated counts and takes the logical union of
absence flags. At the root the stationary-weighted mass of acceptable
buckets is the retention probability. For an observed recipient in a
contact copy, count the donor twice if both tips are observed, and once
if only the recipient is observed. The recipient's inherited state then
contributes no count. Rate and contact probabilities are averaged before
taking the ratio. The tests enumerate retained patterns with both
thresholds, varied masks, and contact components.

## Geographic Gaussian pruning

The Brownian-on-tree likelihood uses independent contrasts: at each internal node
the children's location estimates are combined by precision weighting, the added
variance is the harmonic combination of the child branch variances, and each
child-pair contrast contributes a Gaussian density. This integrates out all
internal locations in closed form and gives the exact likelihood of contrasts under the stated planar
Gaussian model, the continuous counterpart of the discrete
pruning recursion.


# Chapter 6: Sampler correctness

The sampler is a Metropolis-within-Gibbs cycle over scalar parameters and tree
structure. This chapter states the proposal ratios and describes the tests
that compare production moves with exact targets.

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

Each move proposes a uniform increment. For an internal non-root node,
the window is a fraction of the interval from its oldest child to its
parent. Both bounds remain fixed during that move. For the root, the window depends on its oldest child. A dated tip
uses a fixed window from its date density. A proposal that crosses a child
or parent age is rejected. Because the window is unchanged by an accepted
move, the reverse increment has the same density.

## Root scale move

All internal node ages are multiplied by a common factor `exp(lambda(u-1/2))`.
Scaling k independent coordinates contributes a Jacobian of `factor^k`, so the
Hastings term is `k log(factor)`. Invalid results (negative branch) are reverted
and rejected.

## Topology moves

**NNI.** A node's sibling is swapped with one of the node's children. Under
symmetric selection the move is its own inverse, so the Hastings term is zero. The
move is rejected if it would make a child older than its parent.

**Time-preserving SPR.** Choose a prunable child. Detach it and collapse
its parent, keeping the parent's age as a fixed proposal coordinate.
Choose uniformly among edges of the detached backbone crossing that age.
The original sibling edge remains available as a self-transition. In the
reverse direction the same child produces the same detached backbone and
age, so the target-edge sets and selection probabilities agree. The move
preserves the number of prunable children and has no height Jacobian.

The ranked-tree test enumerates all 18 four-tip states at three fixed
event ages. Production NNI and SPR proposals with the exact lexical
target visit all 18 states and meet a fixed frequency tolerance.
This conditional test does not establish mixing at full study scale.

## Checkpointing and exact resumption

The chain's RNG is a serializable xoshiro256** with a four-word state. A
checkpoint writes the RNG words, tree, parameters, frozen step sizes, iteration
counters, and retained draws atomically via a temp-file rename. Resuming
reconstructs an identical chain. A chain interrupted and resumed partway produces
a trajectory bit-identical to an uninterrupted run: same final log posterior,
same topology, identical retained draws, and no duplicated retained draw (gate G12,
`tests/test_recovery.py`). Adaptation of step sizes runs only during warmup and is
frozen before any retained sample, so the retained chain targets the stationary
posterior.

## Coupled rate and time proposals

The rate-clock ridge proposal multiplies clock by f and divides gain and
loss by f, preserving effective lexical rates. Its log Hastings term is
-log(f), from the product of the original coordinates. The independent
conditional-lognormal test checks the induced clock distribution.

The coupled time-clock move multiplies k internal ages by f and divides
clock by f. Its log Hastings term is (k-1) log(f). Dated tips stay fixed,
so effective branch lengths need not be invariant. The acceptance ratio
reevaluates the complete dated target. The dedicated root move scales the
positive gap above the root's oldest child and contributes log(f).

## Composed finite posterior check

`make sampler-check` enumerates 72 conditional states for each required
observation model, combining ranked four-tip topologies, rate settings and
two historical-tip ages. The reference sums complete binary node states
and contact histories, applies the observed-mask retention event, integrates
a dense Gaussian geographic precision matrix and uses independent SciPy
prior and terminal-date densities. It shares no pruning or geographic
message implementation with the production score.

The production acceptance method samples each finite target for 200000
retained steps after 1000 warmup steps. The predeclared maximum frequency
deviation is 0.01 and the log-posterior tolerance is 1e-9. Passing this audit
validates the composed conditional targets. Continuous proposal ratios,
checkpoint recovery and full-data mixing require their separate checks.

## Resampled attachment age

The optional native SPR move prunes the same kind of child as fixed-height
SPR. It suppresses the binary parent. Each backbone edge with enough time
above the detached child is eligible. The proposal chooses an edge uniformly
and samples the reattached parent's age uniformly between the oldest child
and the edge's upper endpoint. The original edge remains eligible.

Detaching the same child after the move yields the same backbone. The number
of prunable children and the eligible edge set agree in reverse. If the old
and new age spans are w and w', the age density ratio is w'/w. The accepted
proposal contributes log(w') minus log(w). A four-tip test fixes the root
and samples uniformly over the valid internal heights. Each balanced topology
has a square of height volume. Each pectinate topology has half a square.
Sampled frequencies meet the declared tolerance for all fifteen labeled
rooted topologies.

## Weighted attachment candidates

An optional attachment move draws four candidates on the detached backbone
and adds the current state to that set. Given the detached child and
backbone, it selects an entry with weight proportional to its full posterior
density divided by its edge and age proposal density. In the augmented
distribution, the selected position has the target density and the other
positions have the proposal density. Resampling the position from those
weights preserves the target marginal.

The parsimony option uses two-state inside and outside messages to compute
the minimum lexical changes at each potential insertion edge. Those counts
use only the fixed columns and the shared backbone and subtree. Missing
states allow either value. The finite-edge probabilities retain a uniform
component and the root stem retains its own mass. Enumerating states on
small trees checks the placement scores. A separate four-tip test checks
the stationary topology frequencies and unbounded root-age distribution
with and without the parsimony weights.

The likelihood guide uses inside and outside lexical messages to score three
age bands per finite edge. Within a site, each insertion score is divided
by the corresponding backbone likelihood. The backbone likelihood and
the normalized detached message introduce only constants common to all
placements. The guide gives each band nonzero probability and uses a
uniform density inside that band. It omits the retention, geography and
contact terms from the guide score. The selected state's importance weight
still uses the complete target and the exact continuous proposal density.

## Joint calendar scale and time block

The calendar-scale move multiplies the inferred ages by exp(t) and divides
the clock and any active diffusion rate by exp(t). Its volume factor is
exp((d-r)t) for d free ages and r reciprocally scaled rates. Random placement
of a unit interval on the log scale and subsequent shrinkage
produce a reversible slice kernel for the calendar, coalescent and rate
priors times that volume factor. The full posterior acceptance calculation
subtracts the prior-density change accounted for by the slice kernel.
[Slice sampling](https://arxiv.org/abs/physics/0009028)

The joint time and clock block uses log gaps above the oldest child for
internal nodes and log gaps above the date offset for sampled tips. The log
clock is a separate coordinate. Computing ages from younger to older gives
a triangular change of variables with log determinant equal to the sum of
those log coordinates. A Gaussian proposal uses the lexical, coalescent and
terminal date gradients to guide its mean. Acceptance uses the full posterior,
both Gaussian densities and the log determinant difference. The step size
adapts during warmup and is saved in the checkpoint.
[Langevin proposal reference](https://doi.org/10.2307/3318418)

With multiple steps the same chart supports a Hamiltonian proposal. Fixed
metric drifts and position-dependent momentum kicks each preserve volume.
The symmetric sequence followed by a momentum flip is its own inverse.
The acceptance rule adds the initial minus final Gaussian kinetic energy to
the full posterior and chart determinant changes. This remains valid when
the deterministic force omits a posterior term. A test checks the inverse
even for a nonconservative force, and a separate run resumes the complete
multi-step kernel exactly.

The lexical gradient combines derivatives of the observed pattern and the
retention probability. The pattern derivative pairs a descendant likelihood
with its complement on each edge. The retention derivative carries the
same saturated count states as the likelihood, along with outside weights
for each count bucket. Rate and contact mixtures average component
derivatives by their mass. Tests compare the result with finite differences
of the production likelihood and check the time-coordinate determinant
and exact resume. The observed pilots measure mixing separately.
[Sampler pilot results](../results/sampler-pilots.md)

## Lexical replica exchange

For each replica with inverse temperature beta, its target is proportional
to the lexical likelihood to the power beta times the other posterior factors.
The temperature-one marginal is the original posterior. The joint target is
the product across replicas. Swapping two complete states leaves the product
of all nonlexical factors fixed. With lexical scores L_a and L_b, the log
exchange ratio is (beta_a - beta_b)(L_b - L_a). Each local transition and each
adjacent exchange preserves the joint target after warmup.

A finite test computes the entire two-replica joint distribution directly and
compares its probabilities with sampled frequencies. The production checkpoint
stores all local RNGs, the exchange RNG, the current walker order, swap counts
and cold draws in one atomic file. An interrupted run resumes exactly at a
nonexchange iteration, including all hot states and cold retained iterations.


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

Protocol 3 generated and fitted 200 datasets with the same serial tree,
binary selection, date and candidate-copy models. Four chains per dataset
retained 1280 root ages. Their nominal 90% intervals covered 177 of 200
true ages. Nontrivial clade recovery averaged 24.7%. The matched run also
fell below the declared coverage target. Its different truth design and
chain budgets limit comparisons with protocol 2. [Recovery evidence](../results/recovery_200.md)

Protocol 4 repeated the same 200 matched datasets with the joint age and
clock move and double the sampling budget. Nominal intervals covered 170
true root ages (85%). Mean clade recovery rose to 28.8% and median root bulk
ESS rose from 20.7 to 176.7. Nine previously covered cases lost coverage
and two gained it. The paired result combines the sampler and budget
changes. The independent replay matched the generator probes and
recomputed the retained interval and clade summaries.
[Protocol 4 evidence](../results/recovery-v4-matched-verification.json)


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


## Full inference budgets

All 16 chains of the protocol 3 inference campaign completed at the declared
budget. Every model failed convergence. Root split R-hat ranged from 1.53 to
1.99 and root bulk ESS from 5.9 to 7.6. The maximum between-chain frequency
difference for a clade was 1.00 in each model. Full checkpoints and the
source-bound replay checks passed. [Full inference results](../results/full-inference.md)

Shorter dated pilots isolated calendar moves and topology transitions at a
common budget. Holding two reference topologies fixed, the four-step
Hamiltonian move increased root bulk ESS from 8.8 and 15.8 to 79.4 and
27.5, respectively. The free-tree pilot improved root split R-hat from 2.115
to 1.336. Adding a lexical likelihood guide brought it to 1.324, with bulk
ESS 9.57 and tail ESS 11.86. The maximum between-chain clade frequency
difference was 1.00 in both free-tree pilots. All runs passed replay, and the
full inference gate remains unmet. [Sampler pilot results](../results/sampler-pilots.md)

The alternative-history registry counts declared age windows and exact
catalog-defined clades over every saved draw. It retains source and context
digests and the counts by chain. The Germanic group occurs in every draw of
all four bundles. The Slavic group occurs in 83.9% of the borrowing draws
and every draw of the other three bundles. These frequencies describe a
specific frozen cohort and unconverged samples. [Observed scenarios](../results/observed-scenarios.json)

The observed prediction experiment withheld 34 concepts or 30 varieties.
The geography fit had a higher saved-draw score on withheld varieties than
tree-only. The comparison across withheld concepts was near zero. The
concept bootstrap does not propagate poor chain mixing. [Held-out results](../results/observed-campaign.md)


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
their full date ranges. The export also contains location marginals for up to
24 recurring clades, conditional on each clade occurring in a retained tree.
Regional fits use inward and outward state messages to draw the root and
clade categories. The workbench places categorical counts at named centroids
and identifies the size of the filtered draw set.

None of the geographic diagnostics establish a homeland. Coordinate spread,
between-chain agreement, and stability to the reference-point model remain
empirical checks.

On known trees the spherical mesh scored better under spherical diffusion
truth. The regional-truth comparison was unresolved. This check evaluates
the spatial processes with tree and time fixed. [Spatial study](../results/spatial-study.md)


# Chapter 12: Sensitivity and robustness

The synthetic sensitivity diagnostic defines fifteen settings before
examining their results. It varies coalescent size, ascertainment, singleton
eligibility, loan filtering, an attestation subset and removal of calibration.
It also compares rate heterogeneity, a relaxed clock, planar geography,
two region grids and two contact priors.

Every configuration keeps its model definition, attestation count and
root quantiles. The half-attestation lever removes source records before
deriving the calibration. The subset can retain the oldest observation, so
a small shift would not imply robustness to different calendar evidence.

The two-chain diagnostics use a known-truth corpus. The subsequent observed
campaign completed all 20 predeclared settings with four full-budget chains
each. All 20 archived fits replayed, and all failed convergence. Eight
settings shifted the saved root mean by more than 15% of baseline. The
changes combine differences in assumptions with poor chain mixing.
[Observed sensitivity results](../results/observed-campaign.md)


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


# Chapter 14: Limitations

The binary transformation discards some dependence among synonyms and
partial judgments. Conditioning on observed masks removes missingness from
the predictive target. Lexical coverage can still be informative about source
sampling and model mismatch. Comparing selected binary patterns also leaves
prediction of wholly new cognate classes outside the held-out language
design.

The contact model considers a finite candidate set and one terminal copy
per character. It can miss documented transfers between distinct classes.
The geographic model conditions on modern reference positions and does
not infer a population boundary. A planar approximation measures degrees
rather than physical distance. The discrete model depends on its chosen
regions and currently uses hard assignments for scored tips.

The dated model uses historical varieties as terminal observations and can
neither infer sampled ancestors nor adjudicate their historical status.
The tree prior and substitution-rate priors help determine time scale.
The independently sourced Latin inscriptions are useful historical
observations with strong geographic and register bias. They do not
calibrate the fitted IE-CoR models.

The observed chains missed convergence, and nominal intervals undercovered
in synthetic recovery. These failures leave historical tree, time and
location comparisons open.
