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
