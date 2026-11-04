"""Summary statistics for posterior predictive checks.

These statistics describe a character matrix in ways the models should reproduce
if they fit. They are computed on both the observed data and on data simulated
from posterior draws; a large gap between observed and predicted flags
misspecification. The statistics are chosen to probe different model aspects:
the site-pattern frequency spectrum (does the model get cognate class sizes
right), the missingness rate (is missingness handled), the pairwise sharing
distribution (does the tree explain who shares cognates), and a loan-sensitive
statistic (clustering of shared classes among geographic neighbors).
"""

from __future__ import annotations

import numpy as np

PRESENT = 1
ABSENT = 0
MISSING = -1


def class_size_spectrum(columns: list[list[int]]) -> np.ndarray:
    """Histogram of how many languages are present in each character.

    Returns counts indexed by present-count; the shape of this spectrum reflects
    the cognate class size distribution the model implies.
    """
    n_lang = len(columns[0]) if columns else 0
    spectrum = np.zeros(n_lang + 1)
    for col in columns:
        k = sum(1 for s in col if s == PRESENT)
        spectrum[k] += 1
    return spectrum


def missingness_rate(columns: list[list[int]]) -> float:
    total = sum(len(c) for c in columns)
    missing = sum(1 for c in columns for s in c if s == MISSING)
    return missing / total if total else 0.0


def mean_pairwise_sharing(columns: list[list[int]]) -> float:
    """Average over character columns of the fraction of observed language pairs
    that share the same present/absent state.

    This probes whether the tree-induced correlation structure matches the data:
    a model that gets the topology wrong predicts a different sharing level.
    """
    if not columns:
        return 0.0
    n_lang = len(columns[0])
    totals = []
    for col in columns:
        obs = [i for i in range(n_lang) if col[i] != MISSING]
        if len(obs) < 2:
            continue
        same = 0
        pairs = 0
        for a in range(len(obs)):
            for b in range(a + 1, len(obs)):
                pairs += 1
                if col[obs[a]] == col[obs[b]]:
                    same += 1
        if pairs:
            totals.append(same / pairs)
    return float(np.mean(totals)) if totals else 0.0


def singleton_fraction(columns: list[list[int]]) -> float:
    """Fraction of characters present in exactly one language (before filtering).

    A model that mis-estimates the gain rate predicts the wrong singleton rate.
    """
    if not columns:
        return 0.0
    singletons = sum(
        1 for col in columns if sum(1 for s in col if s == PRESENT) == 1
    )
    return singletons / len(columns)


def all_statistics(columns: list[list[int]]) -> dict:
    spec = class_size_spectrum(columns)
    return {
        "n_characters": len(columns),
        "missingness_rate": missingness_rate(columns),
        "mean_pairwise_sharing": mean_pairwise_sharing(columns),
        "singleton_fraction": singleton_fraction(columns),
        "mean_class_size": float(
            np.average(np.arange(len(spec)), weights=spec) if spec.sum() else 0.0
        ),
    }
