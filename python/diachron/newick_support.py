"""Helper for matching a taxon set against clade-posterior keys.

The topology diagnostics key clades by the frozenset of their member taxa. A
scenario asks for the support of a specific taxon set, which may or may not be an
exact key; this helper looks it up, returning the posterior probability of that
exact bipartition or zero if it never appears in the posterior.
"""

from __future__ import annotations


def clade_support_for(posteriors: dict, taxa: set[str]) -> float:
    """Posterior probability of the clade defined by exactly `taxa`."""
    target = frozenset(taxa)
    return float(posteriors.get(target, 0.0))
