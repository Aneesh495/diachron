# Experiment plan (frozen before full measurement)

This plan is frozen before the full campaign. A change to it requires a new
version and a new campaign that retains the earlier results.

## Hypotheses and endpoints

The predeclared questions and hypotheses are in `docs/thesis/01_research_questions.md`.
Declared scalar endpoints for convergence are root age, gain, loss, and clock
rate. The decisive comparisons are: dated versus tree-only held-out predictive
density (Q1); borrowing-aware versus dated fit on loan-bearing data (Q2);
likelihood with versus without ascertainment (Q3); root age under full, half, and
no calibration (Q4); strict versus relaxed clock divergence times (Q5); location
posterior under varying geographic priors (Q6); between-chain clade agreement
(Q7); estimate shift under single-observation removal (Q8).

## Workloads

- Primary corpus: 112 varieties, 110 concepts, 1576 variable characters, 140
  dated attestations. Labeled synthetic known-truth.
- Second corpus: 36 varieties, 90 concepts, 60 attestations, younger root and its
  own calibration.
- Recovery: 200 independent synthetic datasets spanning 8-20 varieties, 40-90
  concepts, root ages 3000-8000 years, borrowing probability in {0, 0.15, 0.3},
  missingness in {0, 0.08, 0.2}.

## Independent evaluation units

Recovery uses each synthetic dataset as an independent unit. Held-out prediction
holds out whole concepts (no cognate class spans the partition); masks are applied
before fitting. Chains are the unit for convergence; at least four independently
initialized chains per required model.

## Tolerances and decision rules

- Likelihood oracle (G02): max absolute log-likelihood difference below 1e-9 over
  at least 1000 tiny cases.
- Sampler exactness (G03): finite-target frequency max deviation below 0.01;
  Gamma scale-move mean within 0.15 and variance within 1.2 of analytic.
- Ascertainment (G04): analytic p(variable) within 0.005 of the simulated
  fraction over 400000 simulations.
- Convergence (G06): split-Rhat below 1.01 and bulk and tail ESS above 400 for
  every declared endpoint; otherwise the model is reported not converged.
- Recovery (G05): report mean clade recovery and 90 percent root-age credible
  interval coverage; coverage materially below 0.9 is a reported over-confidence.
- A conclusion is reported stable only if it holds across the relevant sensitivity
  levers; otherwise it is labeled prior- or calibration-dependent.

## Seeds, hardware, and ceilings

- Seeds: chains use fixed base seeds (campaign: 1000 + chain index via a fixed
  multiplier); recovery datasets use 20000 + index. All seeds are recorded in the
  evidence bundle.
- Hardware: the diagnostic and full profiles are declared separately. The CPU
  diagnostic profile runs on a laptop (8 cores, 16 GB). The full acceptance
  profile uses longer chains; it is resumable and checkpointed.
- Ceilings: chains are bounded by explicit warmup and iteration counts; no
  interactive or unbounded work. Downloads and paid resources are not initiated.

## Exclusions

Singleton cognate classes are excluded from the likelihood and counted.
Non-finite proposals are rejected. Poor seeds are not dropped after the fact; a
protocol change requires a new version.

## Profiles

- CPU diagnostic: 4 chains, 400 warmup, 1200 iterations, thin 10. For quick
  feedback; cannot set the acceptance flag.
- Full acceptance: 4 chains, 3000 warmup, 12000 iterations, thin 10, into a
  sealed evidence directory, verified independently.

## Protocol 2

The original synthetic profiles remain available as diagnostics. The observed
diagnostic uses frozen IE-CoR v1.2.1 and ABVD v0.3 input bytes, four independently
initialized chains per model, 400 warmup iterations, 1200 sampling iterations,
and thinning by ten. IE-CoR dates are published terminal-tip distributions.
ABVD has a tree-only baseline with no borrowed calibrations. Reports preserve
failed convergence and identify the executable source separately from later
diagnostic recomputation.

The revised known-truth campaign still uses 200 datasets, the original lexical
regimes and fixed seeds. It uses two chains, 400 warmup iterations, 1200 sampling
iterations and thinning by ten per dataset. The dating evidence cycles through
12, 20 and 32 source observations. Starts use the derived calibration midpoint
and have independently sampled topologies. Root coverage uses the 5th and 95th
sample quantiles and clade recovery excludes the all-tip clade. Per-dataset root
samples and true-clade masses support independent recomputation. Failures stay
in the denominator for completion accounting.

## Protocol 3 sampler pilot

The retained target and observed input bytes stay as in protocol 2. The
sampler pilot compares random starts with starts from bootstrapped lexical
Jaccard distances and average linkage. Those distances initialize topology
and do not enter the posterior. Internal time-slide windows depend on the
interval between a node's oldest child and its parent. Additional proposals
move gain, loss and clock together along constant effective rates, and
couple an internal-age scale move to the inverse clock multiplier. The
Jacobian for each proposal is stated in the sampler note and tested
independently.

Use four chains, 400 warmup iterations, 1200 sampling iterations and thinning
by ten for both initializers on dated IE-CoR. Keep base seeds and the root
start at 8000 years. Select the initializer before the full campaign from
these training diagnostics. All scalar endpoints, topology summaries, move
acceptance counts, runtime and source identity remain visible. The full
campaign retains at least four chains, 3000 warmup iterations and 12000
sampling iterations per posterior model. Pilot diagnostics do not satisfy
the full inference gate.

For protocol 3 full runs, compare nontrivial clades whose mean frequency
across chains is at least 0.1. Require the mean between-chain frequency
standard deviation below 0.01 and maximum frequency difference below
0.1, in addition to the declared scalar criteria and topology movement.
This guards against a large collection of rarely seen splits diluting
between-chain disagreement.

The first protocol 3 pilot showed that a global scale proposal moved the
root rarely. Add a dedicated multiplier on the root's positive gap above
its oldest child, with log Hastings term equal to the log multiplier.
Keep the other proposals, target and diagnostic thresholds unchanged.
Use the distance initializer for the full campaign because its pilot
reached much higher lexical posterior density in every chain. Preserve
the random-start pilot and its nonconvergence. Run a second dated pilot
with the root-gap move before recording the full campaign source.

## Observed held-out campaign

Freeze the observed predictive source before execution. Use the same four
observation models, four independent chains per fit, 3000 warmup iterations,
12000 retained-phase iterations and thinning by ten. Split seed zero holds
out one fifth of complete concept identifiers and one fifth of language
identifiers. Construct characters after removing held-out lexical forms.
Keep external terminal dates and modern reference locations available to all
prediction models that use them. This evaluates lexical prediction conditional
on those external covariates, not prediction of an unknown language's location.

Score each held-out character by the log of the posterior mean probability,
using at most twenty evenly spaced saved draws from each chain. Group paired
score differences by concept and report a 1000-resample concept bootstrap
interval. The interval describes concept variation. Unconverged fits remain
unconverged, regardless of the score interval. Save complete model bundles,
partitions, per-character scores and exact checkpoints. Resume requires the
same source, input digest, partition seed, chain budget and initializer.

## Observed sensitivity protocol 1

Run the twenty configurations declared in `sensitivity/run.py` on the pinned
IE-CoR input, with four chains, 3000 warmup iterations, 12000 sampling
iterations and thinning by ten. Retain each complete bundle, transformation
record, diagnostics and checkpoint. The original fifteen-configuration
synthetic result remains a distinct diagnostic.

The observed half-calibration setting replaces every second terminal density
in sorted taxon order with the same positive-truncated Normal(4000, 4000)
prior in years before 2000 CE. The all-broad setting replaces all terminal
densities with that prior. This does not reuse withheld means and does not
make an ancient tip contemporary. The no-calibration setting is the explicitly
undated tree-only baseline. EDH list size has no effect on these terminal
densities, because EDH assignment to model tips remains unresolved.

The alternative-cognacy setting merges all usable source-proposed links within
a meaning concept. This is an extreme alternative analysis. The source ordinal
scale is recorded but is not interpreted as a probability. Links between
incompatible concepts or missing endpoints stay unused and are counted.

The catalog setting conditions on monophyly of the first source-classification
level in IE-CoR. This is a taxonomic assumption under test. It cannot be used
as independent corroboration of those same groups. Constrained distance starts
honor the groups, and proposals outside the support have zero density. Required
four-model inference keeps no taxonomic constraints.

The observed regional alternatives use fixed Atlantic, Central European,
Caucasus, Iranian, North Indian and South Indian anchors. The three-anchor
alternative uses entries one, three and five of that ordered list. These
anchors are model states and do not represent ancient boundaries. Additional
settings use five-degree tip reference error and diffusion prior log-SD three.
Report predictive and posterior changes alongside chain disagreement. A setting
that fails convergence cannot establish a stable historical conclusion.

## Attachment-age proposal pilot

Compare the protocol 3 move set with the same moves plus a native SPR move
that resamples attachment age. Both fits use dated IE-CoR, independent distance
starts, the same four seeds, 1000 warmup steps, 3000 sampling steps and thinning
by ten. The target and pinned observation bytes are shared. Run both from the
same clean source. Preserve every chain and compare topology movement,
between-chain clade differences, scalar traces and elapsed sampling time.

The additional kernel samples an eligible backbone edge and a uniform date
within its feasible span. The reverse uses the same detached backbone. Its
log Hastings term is the log destination span minus the log original span.
The four-tip continuous-volume test checks all 15 rooted labeled topologies
against their analytically integrated probabilities. This pilot leaves the
full-budget convergence criteria in force.

## Lexical tempering pilot

Use the dated IE-CoR target and the protocol 3 moves. Each independent chain
has eight replicas with fixed inverse temperatures 1, .98, .94, .88, .80,
.65, .47 and .25. Only the lexical log likelihood is multiplied by inverse
temperature. Other factors retain their original strength. Attempt disjoint
adjacent exchanges every twenty local steps, alternating parity. Save draws
from temperature one before any exchange on the same iteration.

Use four independent ensembles, 1000 warmup steps and 3000 sampling steps per
replica, thinning the cold trajectories by ten. Keep the same cold seeds as
the attachment-age pilot control. This fixes the same cold budget and more
compute. Compare scalar traces, clade disagreement, adjacent swap rates,
round trips and total elapsed time. Retain all replica states in checkpoints.
The original convergence thresholds govern any subsequent full run.
