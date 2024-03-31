import numpy as np
import pytest

from diachron import _core


@pytest.mark.parametrize("minimum,columns", [
    (0, [[0, 0, 0, 0], [-1, 0, 1, 1], [0, 1, 0, -1]]),
    (1, [[1, 0, 0, 0], [-1, 0, 1, 1], [0, 1, 0, -1]]),
    (2, [[1, 1, 0, 0], [-1, 0, 1, 1], [1, 1, 0, -1]]),
])
@pytest.mark.parametrize("contacts", [[], [(1, 0, .7), (0, 2, .3), (1, 3, .2), (3, 2, 0.)]])
def test_reverse_branch_score_matches_perturbed_pruning(minimum, columns, contacts):
    # A three-child node followed by a two-child node reuses suffix storage
    # across different arities. Branches need not be in traversal order.
    parent = [5, 0, 0, 5, 0, -1]
    lengths = [.38, .7, .82, 1.1, .54, 0.]
    height = [.7, 0., 0., 0., 0., 1.1]
    taxa = [-1, 0, 1, 2, 3, -1]
    args = (height, taxa, 5, .19, .62, columns, [.2, 1.0, 2.4], minimum, contacts)
    value, gradient = _core.binary_branch_gradient(parent, lengths, *args)
    assert value == pytest.approx(sum(_core.binary_character_scores(parent, lengths, *args)), abs=1e-11)
    reference = []
    step = 1e-5
    for edge in range(5):
        before, after = lengths[:], lengths[:]
        before[edge] -= step
        after[edge] += step
        reference.append((sum(_core.binary_character_scores(parent, after, *args))
                          - sum(_core.binary_character_scores(parent, before, *args))) / (2 * step))
    np.testing.assert_allclose(gradient, [*reference, 0.], atol=1e-8)
