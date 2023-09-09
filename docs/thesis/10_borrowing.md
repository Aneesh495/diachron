# Chapter 10: A bounded contact model

Each lexical character begins as a binary pattern generated on the dated tree.
A categorical variable selects either the inherited pattern, with weight 1, or
one directed contact edge with weight equal to its nonnegative intensity. When
an edge is selected, the donor's latent tip state replaces the recipient's.
The candidate network is fixed by a modern-proximity screen. Directions are
hypotheses with separate intensities, independent exponential priors and equal
prior rates.

For an observed recipient, the copied pattern is scored by observing its
state at the donor and marginalizing the recipient's inherited state. If both
tips were observed and differ, the copy component has probability zero. If
the recipient is missing, copying changes no observed state. For the selection
event, an observed recipient adds one count to the donor and stops contributing
its own inherited count. The likelihood sums over contact and rate categories
before conditioning on retention. The no-contact model follows as all
intensities tend to zero.

The enumeration checks include missing donors and competition between two
edges. They validate the calculation, and the synthetic comparison tests
sensitivity to injected sharing on candidate edges. Independent historical
contact remains a separate question. The model permits only one terminal copy
per character, and a repeated match can arise from shared inheritance or a
copy. The model does not identify contact date, chains of transmission, or
contact outside the candidate set.

IE-CoR's loan classes and ABVD's form loan flags are auditable annotations.
The observed likelihood can be recomputed after removing those forms.
The contact component models visible cognate states and carries no likelihood
for the existence of an annotation.
