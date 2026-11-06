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
