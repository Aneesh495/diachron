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
