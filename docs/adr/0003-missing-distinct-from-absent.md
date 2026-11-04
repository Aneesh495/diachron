# ADR 0003: Represent missing and absent as distinct states

## Status
Accepted.

## Context
A (language, concept) cell with no lexical data is epistemically different from a
cell where a cognate class is genuinely absent. Collapsing them is a common
simplification that biases inference toward shared absences.

## Decision
Use three cell states: present (1), absent (0), and missing (-1). The pruning
likelihood marginalizes a missing tip over all states and pins an absent tip to
0. The binding's negative-state convention carries missing across to Rust
unchanged.

## Consequences
Missing and absent produce distinct, correctly ordered likelihoods (missing is
strictly larger, since it marginalizes), which the tests enforce (gate G04). The
cost is that the transform must know, per concept, which languages have any data,
so it carries a concept-to-languages map; this is cheap and is also what the
evidence trail needs.
