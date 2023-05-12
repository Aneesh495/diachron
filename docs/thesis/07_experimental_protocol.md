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
