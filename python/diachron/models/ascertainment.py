"""Ascertainment correction for variable-character selection.

Diachron only keeps characters that vary among the observed languages: a cognate
class present in every language, or in none, is excluded (see
`io.characters.build_character_matrix`). That selection biases the likelihood,
because the model would otherwise assign some probability to the constant
patterns we never admit. The correction conditions each kept character on the
event that it is variable.

For a single character the conditional likelihood is

    p(pattern | variable) = p(pattern) / p(variable)
                          = p(pattern) / (1 - p(all-absent) - p(all-present)),

so in log space we subtract `log p(variable)` from each character's log
likelihood, where `p(all-absent)` and `p(all-present)` are computed under the
same tree and model by running the pruning likelihood on the constant patterns.

This correction matches the selection process actually used (drop constant
sites). It is validated in `tests/test_ascertainment.py` by simulating the
selection and confirming the corrected likelihood is the right conditional
(gate G04). A generic correction for a different ascertainment scheme would be
wrong here, which is exactly what the brief warns against.
"""

from __future__ import annotations

import math

from diachron import _core


def log_prob_constant_patterns(
    tree,
    model_kind: str,
    model_params: list[float],
    root_prior: list[float],
    n_states: int,
) -> float:
    """Return log p(character is constant) = log(p(all state 0) + p(all 1) + ...).

    For a binary model the constant patterns are all-absent and all-present. For
    a k-state model every single-state constant pattern counts. We evaluate each
    constant pattern with the pruning likelihood and sum the probabilities.
    """
    parent, bl, h, tax, root = tree.to_core_arrays()
    n_tips = sum(1 for t in tax if t >= 0)
    total = 0.0
    for s in range(n_states):
        column = [s] * n_tips
        ll = _core.char_log_likelihood(
            parent, bl, h, tax, root, model_kind, model_params, column, root_prior
        )
        total += math.exp(ll)
    if total <= 0.0:
        return float("-inf")
    return math.log(total)


def ascertainment_log_correction(
    tree,
    model_kind: str,
    model_params: list[float],
    root_prior: list[float],
    n_states: int,
    n_variable_characters: int,
) -> float:
    """Total correction to subtract from the summed character log likelihood.

    Each of the `n_variable_characters` kept characters is conditioned on being
    variable, so we subtract `n * log p(variable)` where
    `p(variable) = 1 - p(constant)`.
    """
    log_p_constant = log_prob_constant_patterns(
        tree, model_kind, model_params, root_prior, n_states
    )
    p_constant = math.exp(log_p_constant) if log_p_constant > float("-inf") else 0.0
    p_variable = 1.0 - p_constant
    if p_variable <= 0.0:
        return float("-inf")
    return -n_variable_characters * math.log(p_variable)
