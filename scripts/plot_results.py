"""Render figures directly from the published diagnostic records."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "docs/results"


def wilson(hits, total):
    z = 1.959963984540054
    fraction = hits / total
    denominator = 1 + z * z / total
    center = (fraction + z * z / (2 * total)) / denominator
    radius = z * np.sqrt(fraction * (1 - fraction) / total + z * z / (4 * total * total)) / denominator
    return center - radius, center + radius


def recovery():
    with gzip.open(RESULTS / "recovery-v2-rows.jsonl.gz", "rt") as source:
        rows = [json.loads(line) for line in source]
    proportions, lower, upper, labels = [], [], [], []
    for value in (0.0, 0.15, 0.3):
        selected = [row for row in rows if row["regime"]["borrow_prob"] == value]
        covered = sum(row["root_in_90ci"] for row in selected)
        proportion = covered / len(selected)
        lo, hi = wilson(covered, len(selected))
        proportions.append(proportion)
        lower.append(proportion - lo)
        upper.append(hi - proportion)
        labels.append(f"{value:g}\n(n={len(selected)})")

    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    ax.errorbar(range(3), proportions, yerr=[lower, upper], fmt="o", capsize=5, color="#245e83")
    ax.axhline(0.9, linestyle="--", color="#885324", label="Nominal interval coverage")
    ax.set(xticks=range(3), xticklabels=labels, ylim=(0.3, 1.0),
           xlabel="Injected borrowing probability", ylabel="Root age interval coverage",
           title="Synthetic recovery across 200 datasets")
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(RESULTS / "recovery-coverage.svg")
    fig.savefig(RESULTS / "recovery-coverage.png", dpi=160)
    plt.close(fig)


def prediction():
    report = json.loads((RESULTS / "predictive-v2.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.4), layout="constrained", sharex=True)
    names = ("dated", "geography", "borrowing")
    for ax, (split, experiment) in zip(axes, report["experiments"].items()):
        means, low, high = [], [], []
        for model in names:
            result = experiment["models"][model]
            mean = result["mean_concept_delta_vs_tree"]
            lo, hi = result["concept_bootstrap_interval_95"]
            means.append(mean)
            low.append(mean - lo)
            high.append(hi - mean)
        ax.errorbar(means, range(3), xerr=[low, high], fmt="o", capsize=4, color="#245e83")
        ax.axvline(0, linestyle="--", color="#777777")
        ax.set(yticks=range(3), yticklabels=names, title=f"Held-out {split}s",
               xlabel="Mean concept log-score delta vs tree baseline")
        ax.invert_yaxis()
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Synthetic predictive diagnostics with concept bootstrap intervals")
    fig.savefig(RESULTS / "predictive-comparison.svg")
    fig.savefig(RESULTS / "predictive-comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    recovery()
    prediction()
