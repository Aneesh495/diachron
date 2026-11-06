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
