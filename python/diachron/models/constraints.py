"""Explicit monophyly constraints for catalog sensitivity experiments."""

from __future__ import annotations


def validate_clades(clades: list[list[int]], n_tips: int) -> list[int]:
    masks = []
    for clade in clades:
        if (
            not 2 <= len(clade) < n_tips or len(set(clade)) != len(clade)
            or any(type(taxon) is not int or not 0 <= taxon < n_tips for taxon in clade)
        ):
            raise ValueError("clade constraints require distinct in-range taxa and nontrivial clades")
        mask = sum(1 << taxon for taxon in clade)
        if mask in masks:
            raise ValueError("duplicate clade constraint")
        if any(mask & other and mask & other not in (mask, other) for other in masks):
            raise ValueError("overlapping clade constraints are incompatible")
        masks.append(mask)
    return masks


def clade_masks(tree) -> set[int]:
    """Return exact descendant sets for internal nodes as taxon bit masks."""
    children = tree.children()
    order = [tree.root]
    for node in order:
        order.extend(children[node])
    descendants = [0] * tree.n_nodes
    masks = set()
    for node in reversed(order):
        if tree.is_tip(node):
            descendants[node] = 1 << tree.taxon[node]
        else:
            mask = 0
            for child in children[node]:
                mask |= descendants[child]
            descendants[node] = mask
            masks.add(mask)
    return masks


def satisfies_clades(tree, masks: list[int]) -> bool:
    return not masks or set(masks) <= clade_masks(tree)
