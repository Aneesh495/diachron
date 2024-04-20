from itertools import product

import pytest

from diachron import _core


@pytest.mark.parametrize("detached", [0, 4])
def test_all_backbone_placement_scores_match_state_enumeration(detached):
    parent = [4, 4, 5, 5, 6, 6, -1]
    taxa = [0, 1, 2, 3, -1, -1, -1]
    columns = [[1, 0, 1, 0], [1, 1, -1, 0], [-1, -1, 1, 1], [0, 0, 0, 0]]
    junction = parent[detached]
    sibling = next(node for node, ancestor in enumerate(parent) if ancestor == junction and node != detached)
    backbone = parent.copy()
    backbone[sibling] = parent[junction]
    reachable = [sibling if junction == 6 else 6]
    for node in reachable:
        reachable.extend(child for child, ancestor in enumerate(backbone)
                         if ancestor == node and child != junction)
    expected = [None] * len(parent)
    for edge in reachable:
        placed = backbone.copy()
        placed[junction] = backbone[edge]
        placed[edge] = junction
        expected[edge] = sum(min(
            sum(states[child] != states[ancestor] for child, ancestor in enumerate(placed) if ancestor >= 0)
            for states in product((0, 1), repeat=len(parent))
            if all(taxon < 0 or column[taxon] < 0 or column[taxon] == states[node]
                   for node, taxon in enumerate(taxa))
        ) for column in columns)
    assert _core.parsimony_placement_costs(parent, taxa, 6, detached, columns) == expected


@pytest.mark.parametrize("detached", [0, 4])
def test_likelihood_placement_grid_matches_full_repruning(detached):
    parent = [4, 4, 5, 5, 6, 6, -1]
    heights = [0., .11, 0., .18, .24, .72, 1.6]
    taxa = [0, 1, 2, 3, -1, -1, -1]
    columns = [[1, 0, 1, 0], [1, 1, -1, 0], [-1, -1, 1, 1], [0, 0, 0, 0]]
    gain, loss = .14, .49
    costs = _core.lexical_placement_likelihoods(parent, heights, taxa, 6, detached, gain, loss, columns)
    junction = parent[detached]
    sibling = next(node for node, ancestor in enumerate(parent) if ancestor == junction and node != detached)
    backbone = parent.copy()
    backbone[sibling] = parent[junction]
    root = sibling if junction == 6 else 6
    assert costs[root] is None
    assert costs[detached] is None
    offsets = []
    for child, values in enumerate(costs):
        if values is None:
            continue
        placed = backbone.copy()
        placed[junction] = backbone[child]
        placed[child] = junction
        lower = max(heights[child], heights[detached])
        upper = heights[backbone[child]]
        for band, value in enumerate(values):
            ages = heights.copy()
            ages[junction] = lower + (upper - lower) * (2 * band + 1) / 6
            lengths = [ages[p] - ages[i] if p >= 0 else 0 for i, p in enumerate(placed)]
            direct = sum(_core.binary_character_scores(placed, lengths, ages, taxa, root, gain, loss, columns, [1.], 0))
            offsets.append(value - direct)
    assert len(offsets) >= 6
    assert offsets == pytest.approx([offsets[0]] * len(offsets), abs=1e-12)
