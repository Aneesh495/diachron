# ADR 0006: A bounded candidate borrowing network, not free reticulation

## Status
Accepted.

## Context
Borrowing creates shared cognates that violate a tree. A fully general reticulation
model (reversible-jump over arbitrary network edges) is attractive but its
target and proposal bookkeeping are hard to validate, and it invites overclaiming
a reconstructed contact history.

## Decision
Model borrowing as intensities on a fixed candidate edge set, with an exponential
prior pulling intensities toward zero. The likelihood lets a recipient match a
donor's class beyond what the tree explains, scored as a mixture.

## Consequences
The model is finite, identifiable, and testable: with no injected borrowing the
intensities stay near zero, and with real shared presence they rise (gate G10).
We avoid reversible-jump bookkeeping we could not independently validate. The
explicit limit, stated in the product and thesis, is that the model can only ask
about contacts on or near its candidate edges; it cannot discover arbitrary
reticulation, and does not claim to.
