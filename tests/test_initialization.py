"""Regression checks for data-derived starts and sampled tip ages."""

import numpy as np

from diachron.infer.initial import distance_starting_tree, lexical_distances


def test_distance_seed_preserves_identical_pairs_and_dated_tips():
    columns = [[1, 1, 0, 0], [0, 0, 1, 1]] * 12
    dates = {0: 3400.0, 2: 1200.0}
    tree = distance_starting_tree(
        ["a", "b", "c", "d"], columns, 8000.0, np.random.default_rng(4), dates,
    )
    assert tree.parent[0] == tree.parent[1]
    assert tree.parent[2] == tree.parent[3]
    assert tree.heights[0] == 3400.0 and tree.heights[2] == 1200.0
    assert tree.root_age() == 8000.0 and tree.times_valid()


def test_unobserved_taxon_distance_is_finite_and_uninformative():
    distance = lexical_distances([[1, 0, -1], [0, 1, -1]])
    assert np.array_equal(distance, [[0, 1, 1], [1, 0, 1], [1, 1, 0]])
