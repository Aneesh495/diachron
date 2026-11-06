"""Conditional simulation follows the complete independently enumerated law."""

from collections import Counter, defaultdict
from itertools import product
import math

import numpy as np
import pytest

from diachron.models.conditional import sample_selected_patterns
from diachron.models.tree import PhyloTree


def test_conditional_copy_and_rate_mixture_matches_complete_enumeration():
    tree = PhyloTree([4, 4, 4, 4, -1], [.2, .4, 0, .6, 1.5], [0, 1, 2, 3, -1])
    gain, loss = .3, .6
    rates = (.5, 1.8)
    edges = ((0, 1, .7), (2, 0, .3))
    pi = (loss / (gain + loss), gain / (gain + loss))
    rng = np.random.default_rng(7419)
    for mask in ([True] * 4, [False, True, True, True]):
        expected = defaultdict(float)
        for rate, inherited, root in product(rates, product((0, 1), repeat=4), (0, 1)):
            mass = pi[root]
            for state, duration in zip(inherited, tree.branch_lengths()):
                retained = math.exp(-(gain + loss) * rate * duration)
                mass *= pi[state] * (1 - retained) + (root == state) * retained
            for donor, recipient, weight in ((0, 0, 1.), *edges):
                visible = list(inherited)
                visible[recipient] = inherited[donor]
                outcome = tuple(state if observed else -1 for state, observed in zip(visible, mask))
                if outcome.count(1) >= 2 and 0 in outcome:
                    expected[outcome] += mass * weight
        n = 30000
        draws = sample_selected_patterns(tree, gain, loss, [mask] * n, 2, rng,
                                         rate_categories=rates, contacts=edges)
        counts = Counter(map(tuple, draws))
        assert set(counts) == set(expected)
        total = sum(expected.values())
        assert max(abs(counts[key] / n - value / total) for key, value in expected.items()) < .012


def test_rare_retention_needs_no_rejection_loop():
    tree = PhyloTree([4, 4, 4, 4, -1], [0, 0, 0, 0, 7000], [0, 1, 2, 3, -1])
    rng = np.random.default_rng(233)
    draws = sample_selected_patterns(tree, 1e-10, .001, [[True] * 4] * 64, 2, rng)
    assert len(draws) == 64
    assert all(draw.count(1) >= 2 and 0 in draw for draw in draws)
    with pytest.raises(ValueError, match="conditional simulation mask"):
        sample_selected_patterns(tree, .0005, .001, [[True, False, False, False]], 2, rng)
