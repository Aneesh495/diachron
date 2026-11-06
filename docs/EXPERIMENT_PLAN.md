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

## Matched binary recovery protocol 3

Keep the protocol 2 recovery result and its coverage deficit. Run a separate
matched experiment to distinguish inference behavior under its own model from
behavior under the lexical simulator and aggregate calibration used earlier.
Use 200 independently seeded datasets, 8 to 19 tips and 40 to 89 retained
characters. Each fit runs four chains with 800 warmup iterations, 3200 sampling
iterations and thinning by ten. Use independent distance starts at 8000 years.
The observed study and its acceptance thresholds retain their own protocol.

Draw dated tip ages from the same positive-truncated normal distributions used
in fitting, with zero, two or four dated tips. Conditional on those ages, simulate
the serial coalescent and draw rates from the declared priors. Cycle through
three observed-mask regimes. Half the design alternates into a finite directed
copy model on a fixed candidate pair. Draw intensities from its fitted prior.
Generate each binary character directly conditional on its mask and retention
rule. The experiment varies tree shape, time scale, copy intensity and missingness.
It does not use source lexical records or synthetic historical claims.

Record every dataset, raw generated patterns and truth, all retained root ages,
true clade frequencies and chain seeds. Recompute 90% interval coverage,
nontrivial clade recovery and the root rank histogram. The rank histogram is
diagnostic because draws are serially dependent. Preserve failures in the
completion count. Report the matched and earlier misspecified designs separately.

## Informed attachment pilot

Keep the dated IE-CoR input, four seeds and 1000/3000/10 budgets from the
attachment-age pilot. Add one kernel to the control move schedule with
`--attachment-resampling`. Each call keeps the current attachment and draws
four independent edge and age alternatives on the detached backbone. The
proposal uses observed-mask lexical distances and positive uniform mixture
components. The choice weights use the full target divided by proposal
density. No proposal parameter is fitted to the pilot results.

Conditional on the selected child and labeled backbone, assign one pool slot
to the target density and draw the other four from the proposal. The
conditional slot probability is proportional to target divided by proposal
density. Resampling this slot is a Gibbs step. Selecting the child uniformly
over the fixed node array preserves the selection probability. The root age
and the detached subtree stay fixed. A boundary state or a prune that would
remove the root leaves the chain unchanged. The ordinary topology and age
moves remain available.

Compare posterior score, root split Rhat and ESS, distinct topologies,
between-chain clade disagreement, move frequency and elapsed time. Preserve
all trajectories. An improvement on a pilot does not change full thresholds.

## Informed root-stem pilot

Repeat the same dated study and 1000/3000/10 budgets with the stem option
enabled. The former root child may be pruned, and the reused attachment
junction can be inserted above the backbone root. Give the stem 18% of
proposal mass when ordinary edges are available, with an exponential age
gap of mean half the coalescent scale. Keep all ordinary edges in the pool.
The same detached backbone and stem density appear in both directions.
The pilot uses a strict clock so changing the root does not require relocating
a fixed edge multiplier. Compare root Rhat, ESS, topology and posterior scores
against the fixed-root informed arm and the control. Keep all frozen gates.

## Spatial observation study

Compare a six-region distance-decayed CTMC and a 96-state local diffusion
mesh on the unit sphere. Both predict the same six directional regions.
Treat the 12-tip tree and its times as known. Use two independent truth
regimes. Simulate regional histories from jump times and spherical histories
by tangent Gaussian steps with per-axis variance at most 0.0025 radians².
The latter approximates Brownian motion on the sphere. Its directional mean
is tested against the first spherical harmonic moment.

Normalize each generator's average initial loss of directional correlation to
0.00016 per year. Fit one multiplier from {0.25, 0.5, 1, 2, 4} for each model
and truth regime by likelihood on 40 training trees. In each truth regime,
hold out three fixed tips per tree on 120 other trees. Select held-out tips
before generating locations. Evaluate mean categorical negative log scores
per independent tree and form a paired 95% bootstrap interval with 1,000
resamples. Report 90% root region set coverage and mean set size. A negative
mesh-minus-region log-score difference favors the mesh. An interval spanning
zero leaves the predictive ranking unresolved. The study conditions on known
trees and gives no evidence of a preferred observed historical geography.

## Conditional-time diagnostic

Use the first retained draw of chains one and four from the full dated IE-CoR
fit at `01459af9ab79bca1c97ac2490624223a294292a2`. For each topology require
all its nontrivial clades and run four chains for 1000/3000/10. Vary the
initial root gap above the oldest daughter, gain and clock across the four
starts. Preserve the same lexical observations and priors. Check that each
trajectory stays on the prescribed topology and replays from the saved
conditional model. Report scalar Rhat, ESS and root means separately for the
two topologies. This diagnostic isolates within-topology time mixing and
keeps the full topology convergence criteria in force for acceptance.

## Coupled calendar-scale pilot

Repeat the fixed-topology diagnostic with the same anchors, initial states,
seeds and budgets. Enable one extra move that scales all inferred calendar
ages together and applies reciprocal changes to the clock and, where active,
the geographic diffusion rate. Evaluate the resulting within-topology Rhat,
bulk and tail ESS, and count the proposals and accepts. Also run the free dated
topology pilot with the same four seeds, distance initialization and 1000/3000/10
budget used in the earlier dated control. Record maximum between-chain clade
frequency difference as well as root diagnostics. These pilots are diagnostics.
The full budgets and scientific thresholds still apply.
