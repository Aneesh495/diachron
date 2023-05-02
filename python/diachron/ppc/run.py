"""Run posterior predictive checks on a saved bundle.

For each of several posterior draws we simulate a fresh character matrix down the
sampled tree with the sampled rates, compute the summary statistics, and compare
the distribution of predicted statistics to the observed value. A Bayesian
p-value near 0 or 1 for a statistic indicates the model cannot reproduce that
feature of the data (gate G07 support, and a direct misspecification probe).

Simulation reuses the forward model from the simulator but conditions on the
exact sampled tree and rates rather than drawing a fresh tree, so the check asks
"could this fitted model have produced data like ours".
"""

from __future__ import annotations

import numpy as np

from ..infer.artifacts import PosteriorBundle
from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study
from ..models.tree import PhyloTree
from . import stats


def _parse_newick_to_tree(newick: str, tip_labels: list[str]) -> PhyloTree:
    """Reconstruct a PhyloTree (structure + branch lengths) from a Newick string.

    We parse the nested structure, assigning taxon indices by matching tip labels
    to their position in `tip_labels`, and recover node heights from branch
    lengths by accumulating from the tips upward.
    """
    label_to_taxon = {lab: i for i, lab in enumerate(tip_labels)}
    nodes_parent: list[int] = []
    nodes_taxon: list[int] = []
    nodes_branch: list[float] = []
    pos = 0

    def new_node(taxon: int, branch: float) -> int:
        nodes_parent.append(-1)
        nodes_taxon.append(taxon)
        nodes_branch.append(branch)
        return len(nodes_parent) - 1

    import re

    def parse() -> int:
        nonlocal pos
        if newick[pos] == "(":
            pos += 1
            children = []
            while True:
                children.append(parse())
                if newick[pos] == ",":
                    pos += 1
                    continue
                if newick[pos] == ")":
                    pos += 1
                    break
            branch = read_length()
            node = new_node(-1, branch)
            for c in children:
                nodes_parent[c] = node
            return node
        else:
            m = re.match(r"[^():,;]+", newick[pos:])
            token = m.group(0)
            pos += len(token)
            label = token.split(":")[0]
            branch = 0.0
            if ":" in token:
                branch = float(token.split(":")[1])
            elif pos < len(newick) and newick[pos] == ":":
                branch = read_length()
            return new_node(label_to_taxon.get(label, -1), branch)

    def read_length() -> float:
        nonlocal pos
        if pos < len(newick) and newick[pos] == ":":
            m = re.match(r":[-0-9.eE+]+", newick[pos:])
            if m:
                pos += len(m.group(0))
                return float(m.group(0)[1:])
        return 0.0

    root = parse()

    n = len(nodes_parent)
    heights = [0.0] * n
    # Height of a node = max over children of (child height + child branch).
    children: dict[int, list[int]] = {i: [] for i in range(n)}
    for i, p in enumerate(nodes_parent):
        if p >= 0:
            children[p].append(i)

    order: list[int] = []
    stack = [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            order.append(node)
        else:
            stack.append((node, True))
            for c in children[node]:
                stack.append((c, False))
    for node in order:
        if children[node]:
            heights[node] = max(heights[c] + nodes_branch[c] for c in children[node])

    labels = {i: tip_labels[nodes_taxon[i]] for i in range(n) if nodes_taxon[i] >= 0}
    return PhyloTree(parent=nodes_parent, heights=heights, taxon=nodes_taxon, tip_labels=labels)


def _simulate_binary_columns(tree: PhyloTree, gain: float, loss: float,
                             n_chars: int, rng: np.random.Generator) -> list[list[int]]:
    """Forward-simulate binary characters down the tree under the gain/loss model."""
    z = gain + loss
    pi_present = gain / z
    ch = tree.children()
    bl = tree.branch_lengths()
    root = tree.root
    n_tips = sum(1 for i in range(tree.n_nodes) if tree.is_tip(i))
    columns = []
    for _ in range(n_chars):
        state = {root: 1 if rng.random() < pi_present else 0}
        stack = list(ch[root])
        while stack:
            node = stack.pop()
            p = tree.parent[node]
            t = bl[node]
            e = np.exp(-z * t)
            if state[p] == 0:
                p_present = pi_present * (1 - e)
            else:
                p_present = pi_present + (loss / z) * e
            state[node] = 1 if rng.random() < p_present else 0
            stack.extend(ch[node])
        col = [-1] * n_tips
        for i in range(tree.n_nodes):
            if tree.is_tip(i):
                col[tree.taxon[i]] = state[i]
        columns.append(col)
    return columns


def run_posterior_predictive(bundle_dir: str, n_draws: int = 40, seed: int = 0) -> dict:
    bundle = PosteriorBundle.load(bundle_dir)
    corpus = build_primary_corpus()
    study = prepare_study(corpus.dataset, corpus.attestations)
    observed = stats.all_statistics(study.columns)

    rng = np.random.default_rng(seed)
    # Pool draws across chains and subsample.
    all_draws = [d for chain in bundle.chains for d in chain]
    if not all_draws:
        return {"error": "empty bundle"}
    picks = rng.choice(len(all_draws), size=min(n_draws, len(all_draws)), replace=False)

    predicted: dict[str, list[float]] = {
        "missingness_rate": [],
        "mean_pairwise_sharing": [],
        "singleton_fraction": [],
        "mean_class_size": [],
    }
    n_chars = len(study.columns)
    for idx in picks:
        d = all_draws[idx]
        tree = _parse_newick_to_tree(d["newick"], bundle.tip_labels)
        cols = _simulate_binary_columns(tree, d["gain"], d["loss"], n_chars, rng)
        s = stats.all_statistics(cols)
        for k in predicted:
            predicted[k].append(s[k])

    # Bayesian p-values: fraction of predicted >= observed.
    pvalues = {}
    for k, vals in predicted.items():
        arr = np.array(vals)
        pvalues[k] = {
            "observed": observed[k],
            "predicted_mean": float(arr.mean()),
            "predicted_sd": float(arr.std()),
            "bayesian_p": float(np.mean(arr >= observed[k])),
        }
    return {
        "gate": "G07-support",
        "model": bundle.model_name,
        "n_draws_used": int(len(picks)),
        "statistics": pvalues,
    }
