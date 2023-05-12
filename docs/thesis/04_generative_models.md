# Chapter 4: Generative models

This chapter writes the complete generative model and its posterior, and shows
the limiting cases in which the full model reduces to a simpler one.

## Unit of inheritance

An observation attaches to a (language, concept) cell. The transform in
`diachron.io.characters` maps cognate judgments to binary presence/absence
characters: one character per cognate class within a concept. For language L and
class K the state is present (L has a lexeme judged in K), absent (L has lexical
data for the concept but no form in K), or missing (L has no lexical data for the
concept). Singletons (present in one language) are excluded as carrying no
internal-structure signal and counted in the ascertainment accounting. Loans are
annotated, not folded silently into inheritance; partial cognates are their own
classes.

## Character evolution

Each class evolves along the tree as a two-state continuous-time Markov chain with
instantaneous rate matrix

    Q = [ -gain   gain ]
        [  loss  -loss ]

so a lineage gains the class at rate `gain` and loses it at rate `loss`. The
stationary distribution is proportional to `(loss, gain)`. The transition matrix
over a branch of length t has the closed form

    P(t) = [ pi0 + pi1 e   pi1 - pi1 e ]
           [ pi0 - pi0 e   pi1 + pi0 e ],   e = exp(-(gain+loss) t),

with `pi0 = loss/(gain+loss)`, `pi1 = gain/(gain+loss)`. This is validated
against the matrix exponential of Q in the Rust test suite.

## Pruning likelihood

For one character the probability of the observed tip states, marginalized over
internal node states, is the Felsenstein pruning likelihood. Define the
conditional likelihood vector `L_v(s) = P(data below v | v in state s)`. Tips are
indicators of the observed state, or the all-ones vector when missing (which is
why missing and absent differ). At an internal node v with children c,

    L_v(s) = prod_c ( sum_j P_sj(t_c) L_c(j) ).

The character log likelihood is `log sum_s pi(s) L_root(s)`. Chapter 5 gives the
log-space scaling that keeps this stable on deep trees and the ascertainment
conditioning.

## Clock and calibration

A strict clock scales branch durations in years to expected substitutions by a
single rate; a relaxed uncorrelated-lognormal clock gives each branch its own
rate from a lognormal with shared mean and log-standard-deviation sigma. As sigma
goes to zero the relaxed clock reduces exactly to the strict clock. Calibrations
place densities on node ages: an aggregate root-age interval derived from the
oldest reliable attestation, which is the minimal structure that makes the strict
clock's time scale identifiable.

## Tree prior

Node times follow a constant-size coalescent: with k lineages the pairwise
coalescent rate is `k(k-1)/(2N)` and the waiting time to the next coalescence is
exponential with that rate. The density is written out in
`diachron.models.priors.coalescent_log_prior`.

## Geography

Lineage location is a latent state. Two models are compared rather than assumed. A
region CTMC moves a lineage among a small set of named regions. A Brownian-on-tree
model diffuses latitude and longitude with variance sigma^2 per year; its marginal
likelihood of the tip coordinates, integrating out internal locations, is the
Gaussian independent-contrasts pruning, the continuous analogue of the character
recursion. Tip coordinates enter as uncertain observations with their own error;
they are modern reference points, never historical boundaries.

## Borrowing

A bounded contact model places a borrowing intensity on each edge of a fixed
candidate network, with an exponential prior pulling intensities toward zero so
the model prefers no borrowing unless the data support it. This is a finite,
identifiable alternative to unrestricted reticulation, which the platform
deliberately does not claim.

## The full posterior

    log p(tree, times, rates, geo, borrow | data)
      = sum_characters log p(character | tree, times, rates)  [pruning]
        - n_variable * log p(variable | model)                 [ascertainment]
      + log p(node ages | calibrations)                        [calibration]
      + log p(tip coords | tree, times, sigma^2)               [geography]
      + log p(borrow intensities)                              [borrowing prior]
      + log p(tree)                                            [coalescent]
      + log p(rates)                                           [clock prior]
      + const.

## Limiting cases

- Borrowing intensities to zero: the borrowing model reduces to the dated model.
- Geographic sigma to infinity: locations become uninformative and the geography
  term is constant, recovering the dated model.
- Relaxed-clock sigma to zero: the relaxed clock reduces to the strict clock.
- No calibration: the calibration term drops and the time scale is identified only
  up to the clock-rate/time trade-off, which Chapter 12 shows empirically.
