"""Diachron: posterior inference over language family history.

The package is organized around a few contracts:

* ``diachron.io`` reads CLDF lexical data and dated attestations into typed,
  validated tables, and transforms cognate judgments into the character matrix
  the likelihood consumes.
* ``diachron.models`` defines the generative models (character evolution,
  clock, tree prior, geography, borrowing) and the composed log posterior.
* ``diachron.infer`` drives the MCMC, checkpointing, and resumption.
* ``diachron.diagnostics`` computes convergence and mixing summaries.
* ``diachron.ppc`` and ``diachron.sensitivity`` implement posterior predictive
  and robustness checks.
* ``diachron.registry`` records alternative historical scenarios against the
  exact model, assumptions, and posterior support that produced them.

The compiled inner loops (pruning likelihood, MCMC moves, geographic kernels)
live in the Rust extension ``diachron_core``.
"""

__version__ = "0.1.0"
