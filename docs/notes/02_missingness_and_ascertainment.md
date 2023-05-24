# Note 2: Missingness and ascertainment

Two selection effects shape the likelihood, and conflating them would bias
inference.

The first is missingness. A missing cell is not an absence: the pruning likelihood
marginalizes a missing tip over all states, while an absent tip is pinned to state
0. The test suite confirms the likelihood of a pattern with a missing tip is
strictly greater than the same pattern with that tip set absent, because
marginalizing admits both possibilities. If missing were silently treated as
absent, the model would overstate how much a sparsely-attested language is known
to lack, and would pull the tree toward spurious shared absences.

The second is ascertainment. Only variable characters are kept, so the kept set is
conditioned on not being constant. The correction divides each character's
likelihood by the probability that a character is variable under the current tree
and model, computed by running the pruning likelihood on the all-absent and
all-present patterns and summing. The validation is a direct simulation: the
analytic variable-probability matches the fraction of forward simulations that
come out variable, to within 0.005 over 400000 draws. This is the specific point
the brief warns about: a generic ascertainment correction for a different
selection scheme would not reproduce the simulated fraction, so matching the
actual selection is what makes the correction correct rather than merely present.

A negative result worth keeping: on the study corpus, turning ascertainment off
did not move the root-age estimate materially (Note 10). That does not make the
correction optional; it changes likelihood values and the relative fit of
constant-heavy versus variable-heavy regimes, which matters for model comparison
even where it does not move a particular scalar endpoint.
