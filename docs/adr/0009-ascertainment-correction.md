# ADR 0009: Condition kept characters on being variable

## Status
Accepted.

## Context
Only variable cognate-class characters are kept (singletons and constant patterns
are excluded). The likelihood must account for this selection, or it will be
biased toward the constant patterns it never admits.

## Decision
Condition each kept character on being variable by subtracting
`log p(variable) = log(1 - p(constant))` from its log likelihood, with
`p(constant)` computed by running the pruning likelihood on the constant patterns
under the current tree and model.

## Consequences
The correction matches the exact selection applied, validated by simulating the
selection and matching the analytic variable-probability to the empirical fraction
within 0.005 (gate G04). The cost is two extra constant-pattern likelihood
evaluations per posterior call; this is bounded and small relative to the full
character block. A generic correction for a different ascertainment scheme would
be wrong, so the correction is tied to the specific drop-constant-sites rule.
