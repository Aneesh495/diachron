"""Pure-Python independent likelihood and posterior oracles for the tests.

These are deliberately simple and slow reference implementations written without
reference to the Rust pruning code. The character likelihood here enumerates all
internal-node state assignments; the finite-target posterior here enumerates a
small discrete space exactly. Agreement between these oracles and the production
Rust/Python code is the evidence for gates G02 and G03.

Keeping the oracle in Python with explicit state enumeration makes the
differential test meaningful: a bug shared between the implementation and its
check would have to be independently reproduced here, which a transcription of
the same algorithm would not guarantee but an independent enumeration does.
"""

from __future__ import annotations

import itertools
import math


def transition_binary(gain: float, loss: float, t: float) -> list[list[float]]:
    """Closed-form 2-state transition matrix, derived independently here."""
    z = gain + loss
    e = math.exp(-z * t)
    pi0 = loss / z
    pi1 = gain / z
    return [
        [pi0 + pi1 * e, pi1 - pi1 * e],
        [pi0 - pi0 * e, pi1 + pi0 * e],
    ]


def enumerate_log_likelihood(
    parent: list[int],
    branch_len: list[float],
    taxon: list[int],
    root: int,
    gain: float,
    loss: float,
    column: list[int],
    root_prior: list[float],
    n_states: int = 2,
) -> float:
    """Likelihood by summing over every internal-node and missing-tip state.

    `column[t]` is the observed state of taxon t, or negative for missing. This
    mirrors nothing in the Rust code; it builds the joint probability directly.
    """
    n_nodes = len(parent)
    trans = [transition_binary(gain, loss, branch_len[i]) for i in range(n_nodes)]

    internal = [i for i in range(n_nodes) if taxon[i] < 0]
    missing_tips = [i for i in range(n_nodes) if taxon[i] >= 0 and column[taxon[i]] < 0]
    free = internal + missing_tips

    total = 0.0
    for assignment in itertools.product(range(n_states), repeat=len(free)):
        state = [-1] * n_nodes
        for idx, node in enumerate(free):
            state[node] = assignment[idx]
        for i in range(n_nodes):
            if taxon[i] >= 0 and column[taxon[i]] >= 0:
                state[i] = column[taxon[i]]
        term = root_prior[state[root]]
        for i in range(n_nodes):
            p = parent[i]
            if p >= 0:
                term *= trans[i][state[p]][state[i]]
        total += term
    return math.log(total) if total > 0 else float("-inf")
