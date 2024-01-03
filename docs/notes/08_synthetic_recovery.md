# Note 8: Synthetic recovery and coverage

The protocol 2 recovery campaign completed 200 known-truth datasets at
`8c47e96bab0f326d8476a90c418304018397df06`. Regimes vary tree size, concept
count, root age, injected borrowing, missingness and synthetic attestation
coverage. Each fit used two chains with 400 warmup iterations and 1200
sampling iterations. The all-tip clade is excluded from topology scoring.

Nominal 90% root intervals covered 144 of 200 true ages, or 72%. Mean
nontrivial clade recovery was 0.4132 at the declared frequency threshold of
0.5. The [summary](../results/recovery-v2.json) and compressed per-dataset
rows preserve every dataset, including its regime, root samples and
true-clade frequencies. These support independent recomputation.

Coverage is below the nominal target. This is a measured failure of interval
calibration at the executed budget. It does not establish which mechanism is
responsible. The simulator, likelihood assumptions, synthetic calibration
construction and short-chain topology mixing can each contribute. Recovery
under a correctly matched simulator and recovery under injected violations
must remain distinguishable.

A posterior interval can become narrow under a misspecified model and still
be wrong. More retained samples can reduce Monte Carlo error while preserving
that bias. Better mixing and a separate matched-model calibration study are
needed before interpreting a coverage deficit as a purely statistical model
failure or a purely computational failure.
