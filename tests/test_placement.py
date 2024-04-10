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
