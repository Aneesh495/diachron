"""Gate G04: the ascertainment correction matches the selection process.

The correction conditions each kept character on being variable (not constant).
We validate it two ways:

1. **Simulated selection.** Simulate characters under the model, keep only the
   variable ones, and estimate the empirical fraction of simulations that are
   variable. The analytic `p(variable)` used by the correction must match that
   empirical fraction. If they disagree, the correction does not describe the
   selection actually applied.

2. **Missing versus absent.** The likelihood of a pattern with a missing tip must
   differ from the same pattern with that tip pinned to absent, and the missing
   version must be larger (it marginalizes rather than committing to a state).
   This is the distinct-likelihood requirement in the gate.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from diachron import _core
from diachron.models.ascertainment import log_prob_constant_patterns
from diachron.models.tree import PhyloTree


def _four_tip_tree():
    # Balanced 4-tip tree with known branch lengths.
    parent = [4, 4, 5, 5, 6, 6, -1]
    heights = [0.0, 0.0, 0.0, 0.0, 0.4, 0.5, 0.9]
    taxon = [0, 1, 2, 3, -1, -1, -1]
    return PhyloTree(parent=parent, heights=heights, taxon=taxon)


def _simulate_pattern(tree, gain, loss, rng):
    """Draw tip states by forward simulation down the tree."""
    z = gain + loss
    pi = [loss / z, gain / z]
    root = tree.root
    bl = tree.branch_lengths()
    state = {root: 0 if rng.random() < pi[0] else 1}
    # Preorder.
    ch = tree.children()
    stack = list(ch[root])
    while stack:
        node = stack.pop()
        p = tree.parent[node]
        t = bl[node]
        e = math.exp(-z * t)
        # P(end present | start) from the closed-form transition.
        if state[p] == 0:
            p_present = (gain / z) * (1 - e)
        else:
            p_present = (gain / z) + (loss / z) * e
        state[node] = 1 if rng.random() < p_present else 0
        stack.extend(ch[node])
    return [state[i] for i in range(tree.n_nodes) if tree.is_tip(i)]


def test_p_variable_matches_simulation():
    tree = _four_tip_tree()
    gain, loss = 0.7, 1.1
    z = gain + loss
    root_prior = [loss / z, gain / z]

    log_p_const = log_prob_constant_patterns(tree, "binary", [gain, loss], root_prior, 2)
    p_const_analytic = math.exp(log_p_const)
    p_var_analytic = 1.0 - p_const_analytic

    rng = np.random.default_rng(13)
    n = 400_000
    variable = 0
    for _ in range(n):
        pat = _simulate_pattern(tree, gain, loss, rng)
        if len(set(pat)) > 1:
            variable += 1
    p_var_empirical = variable / n
    assert abs(p_var_analytic - p_var_empirical) < 0.005, (
        f"analytic {p_var_analytic} vs empirical {p_var_empirical}"
    )


def test_missing_differs_from_absent():
    tree = _four_tip_tree()
    gain, loss = 0.7, 1.1
    z = gain + loss
    root_prior = [loss / z, gain / z]
    parent, bl, h, tax, root = tree.to_core_arrays()

    absent = [1, 0, 1, 0]
    missing = [1, -1, 1, 0]
    ll_absent = _core.char_log_likelihood(
        parent, bl, h, tax, root, "binary", [gain, loss], absent, root_prior
    )
    ll_missing = _core.char_log_likelihood(
        parent, bl, h, tax, root, "binary", [gain, loss], missing, root_prior
    )
    assert ll_missing > ll_absent + 1e-6, f"missing {ll_missing} not > absent {ll_absent}"


def test_changing_states_changes_likelihood():
    """Patterns that are not related by a tree symmetry have distinct likelihoods.

    The balanced tree makes sibling tips exchangeable, so [1,0,1,0] and [1,0,0,1]
    (which swap the two tips under one clade) are genuinely equal: that is a
    correct symmetry, not a collision. We instead compare patterns that break the
    symmetry: grouping cognate presence with one clade versus splitting it across
    clades changes the likelihood.
    """
    tree = _four_tip_tree()
    gain, loss = 0.7, 1.1
    z = gain + loss
    root_prior = [loss / z, gain / z]
    parent, bl, h, tax, root = tree.to_core_arrays()
    # [1,1,0,0] matches clade {0,1}; [1,0,1,0] crosses clades; [1,1,1,0] is three
    # present. These are not symmetry-equivalent and must differ.
    patterns = [[1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 1, 0]]
    lls = [
        _core.char_log_likelihood(
            parent, bl, h, tax, root, "binary", [gain, loss], p, root_prior
        )
        for p in patterns
    ]
    assert len(set(round(x, 9) for x in lls)) == len(lls), f"likelihoods collided: {lls}"

    # And the clade-matching pattern [1,1,0,0] beats the clade-crossing [1,0,1,0]:
    # the tree favors cognate sharing within a clade.
    assert lls[0] > lls[1], f"clade-matching {lls[0]} not > clade-crossing {lls[1]}"


def test_mixture_ascertainment_matches_independently_enumerated_selected_patterns():
    from itertools import product


    parent = [-1, 0, 0, 0, 0]
    lengths = [0.0, 0.3, 0.8, 1.1, 0.9]
    heights = [2.0, 1.7, 1.2, 0.9, 1.1]
    taxon = [-1, 0, 1, 2, 3]
    gain, loss = 0.4, 0.9
    prior = [loss / (gain + loss), gain / (gain + loss)]
    rates = [0.15, 1.85]
    for observed in ((0, 1, 2), (0, 2, 3), (0, 1, 2, 3)):
        for minimum_present in (1, 2):
            columns = []
            probabilities = []
            for pattern in product((0, 1), repeat=len(observed)):
                if sum(pattern) < minimum_present or 0 not in pattern:
                    continue
                column = [-1] * 4
                for tip, state in zip(observed, pattern):
                    column[tip] = state
                columns.append(column)
                # For this star the root is the only unknown ancestral state.
                # Independently sum it out using the closed-form binary CTMC.
                components = []
                for rate in rates:
                    prob = 0.0
                    for ancestor, root_mass in enumerate(prior):
                        mass = root_mass
                        for tip, state in zip(observed, pattern):
                            persist = math.exp(-(gain + loss) * rate * lengths[tip + 1])
                            mass *= prior[state] * (1 - persist) + (ancestor == state) * persist
                        prob += mass
                    components.append(prob)
                probabilities.append(sum(components) / len(components))
            expected = [math.log(p / sum(probabilities)) for p in probabilities]
            actual = _core.binary_character_scores(
                parent, lengths, heights, taxon, 0, gain, loss, columns, rates, minimum_present,
            )
            assert actual == pytest.approx(expected, abs=1e-12)


def test_tree_with_disconnected_cycle_rejected_before_pruning():
    with pytest.raises(ValueError, match="unreachable nodes or a cycle"):
        _core.binary_character_scores(
            [-1, 0, 3, 2], [0.0, 1.0, 1.0, 1.0], [1.0, 0.0, 1.0, 0.0],
            [-1, 0, -1, -1], 0, 0.4, 0.8, [[1]], [1.0], 0,
        )
