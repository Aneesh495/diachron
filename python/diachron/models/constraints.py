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


def satisfies_clades(tree, masks: list[int]) -> bool:
    if not masks:
        return True
    children = tree.children()
    descendants = {}
    pending = [(tree.root, False)]
    while pending:
        node, visited = pending.pop()
        if tree.is_tip(node):
            descendants[node] = 1 << tree.taxon[node]
        elif visited:
            mask = 0
            for child in children[node]:
                mask |= descendants[child]
            descendants[node] = mask
        else:
            pending.append((node, True))
            pending.extend((child, False) for child in children[node])
    return set(masks) <= set(descendants.values())
