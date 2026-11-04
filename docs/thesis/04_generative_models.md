# Chapter 4: Generative models

## Character and rate model

The unit is the membership of a language variety in a cognate class for a
particular concept. One source judgment can make the character present.
Absence requires full classification of the recorded forms in the cell.
A missing lexical cell or a cell with inadequate classification contributes
an unknown state. Loan judgments and multiple or partial classifications
remain in the evidence tables. They do not define mutually exclusive
states for a single meaning.

Let g and l be the gain and loss rates, and c the positive clock
multiplier. The two-state generator has off-diagonals g and l and
stationary distribution (l, g)/(g+l). Branch v has elapsed time
d_v = h_parent - h_v. The transition is exp(Q c d_v r_v a), where r_v
is one for the strict clock and a is one for the homogeneous site model.
Under the relaxed clock, log(r_v) follows a Normal with mean -s²/2 and
standard deviation s for fixed model dispersion s. Under the site
mixture, a is integrated over equal-probability discrete-gamma categories
scaled to arithmetic mean one. The clock multiplier and g and l enter
the likelihood as products. Proper priors regularize that redundancy.

## Observation and posterior

For character j, let x_j be its observed pattern and m_j its
observed-tip mask. Let A_j be the event of at least one observed absence
and either one or two observed presences, depending on the singleton
policy. For a fixed tree, define P_j as the mixture over rate categories
and possible contact edges. A contact component copies one donor's
latent terminal state to the recipient. It is averaged with the inherited
component before dividing by selection probability.

    log target = sum_j [log P_j(x_j | m_j) - log P_j(A_j | m_j)]
               + log p(ordered labeled tree | observed tip ages, N)
               + sum of dated tip and node calibration log densities
               + geographic log contrast or region likelihood
               + log p(g, l, c, r, geographic rate, contact intensities).

The denominator is omitted in the declared uncorrected sensitivity run.
Lognormal priors on g and l have median 0.001 and log standard deviation
1.5. The clock multiplier has median one and log standard deviation 0.5.
The planar geographic diffusion variance has median 0.004 squared degrees
per year and log standard deviation two. Contact weights have independent
exponential priors, with default rate ten. Each prior applies in the model
where its parameter is sampled.

For each waiting interval with k extant backward lineages, the constant-size
coalescent contributes survival under total hazard k(k-1)/(2N). Each labeled
merge contributes its own hazard 1/N. Historical tip dates introduce lineage
sampling events. The IE-CoR model uses terminal-tip densities in years before
2000 CE. The synthetic model can use an aggregate root-age calibration.
The ABVD tree-only analysis does not assign independent historical dates.

The dated baseline omits geography and contact. The geographic variant adds
planar contrasts or a symmetric region CTMC. The borrowing variant includes
contact copying. The contact likelihood approaches the dated lexical
likelihood as intensities approach zero. The relaxed clock approaches the
strict clock as its dispersion approaches zero.
