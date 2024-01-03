# Note 7: Convergence

The declared scalar endpoints are root age, gain, loss and clock rate. Each
must have split and rank-normalized R-hat below 1.01 and bulk and tail ESS
above 400. Independent tests compare the estimators with independent draws,
autoregressive series and deliberately separated chains.

Topology has a separate acceptance condition. Protocol 3 compares nontrivial
clades with pooled sampling frequency at least 0.1. Mean between-chain
frequency standard deviation must be below 0.01 and the maximum chain
frequency difference below 0.1. Every chain must visit multiple topologies.
A count of distinct trees does not prove agreement on high-mass clades.

The original observed diagnostic missed the scalar criteria. The protocol 3
pilots retained all four chains for both random and distance initializers.
Distance starts reached much higher lexical posterior density, but independent
chains still occupied different clade modes. A dedicated root-gap proposal
increased exploration of the root date without resolving that disagreement.
The pilots therefore did not satisfy the inference gate.

The full campaign fixes four chains, 3000 warmup iterations, 12000 sampling
iterations and thinning by ten. The source-bound campaign verifier recomputes
full-precision diagnostics from saved draws, checks complete final checkpoints
and requires unique chain seeds. It preserves failed convergence. A pilot
comparison, an accepted local move and a completed run are distinct forms of
evidence, and none alone establishes stationary sampling.
