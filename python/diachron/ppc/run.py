"""Posterior predictive checks on the exact study and draws saved in a bundle."""

from __future__ import annotations

import numpy as np

from ..infer.artifacts import PosteriorBundle
from diachron import _core
from ..models.tree import PhyloTree
from ..models.conditional import sample_selected_patterns
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


def _replicate(
    tree, gain, loss, rates, observed, minimum, rng, contacts=(), intensities=(), clock=(),
):
    """Rejection sample by observed mask and the same rule used for fitting."""
    n_characters, n_tips = observed.shape
    remaining = np.arange(n_characters)
    replicated = np.full(observed.shape, -1, dtype=np.int8)
    root = tree.root
    children = tree.children()
    order = [root]
    for node in order:
        order.extend(children[node])
    branches = tree.branch_lengths()
    if clock:
        branches = [length * multiplier for length, multiplier in zip(branches, clock)]
    pi = gain / (gain + loss)
    for _ in range(512):
        if not remaining.size:
            return replicated.tolist()
        categories = rng.choice(rates, size=remaining.size)
        values = np.empty((tree.n_nodes, remaining.size), dtype=np.int8)
        values[root] = rng.random(remaining.size) < pi
        leaves = np.empty((remaining.size, n_tips), dtype=np.int8)
        for node in order:
            parent = tree.parent[node]
            if parent >= 0:
                retain = np.exp(-(gain + loss) * categories * branches[node])
                probability = pi * (1.0 - retain) + values[parent] * retain
                values[node] = rng.random(remaining.size) < probability
            if tree.is_tip(node):
                leaves[:, tree.taxon[node]] = values[node]
        if contacts:
            weights = np.array([1.0, *intensities])
            active = rng.choice(len(weights), size=remaining.size, p=weights / weights.sum())
            for edge, (donor, recipient) in enumerate(contacts, start=1):
                sites = active == edge
                leaves[sites, recipient] = leaves[sites, donor]
        mask = observed[remaining]
        if minimum:
            present = ((leaves == 1) & mask).sum(axis=1)
            absent = ((leaves == 0) & mask).any(axis=1)
            accepted = (present >= minimum) & absent
        else:
            accepted = np.ones(remaining.size, dtype=bool)
        replicated[remaining[accepted]] = np.where(mask[accepted], leaves[accepted], -1)
        remaining = remaining[~accepted]
    if remaining.size:
        edges = [(donor, recipient, rate) for (donor, recipient), rate in zip(contacts, intensities)]
        replicated[remaining] = sample_selected_patterns(
            tree, gain, loss, observed[remaining], minimum, rng, rate_categories=rates,
            contacts=edges, clock_multipliers=clock,
        )
    return replicated.tolist()


def _spatial_dispersion(points: np.ndarray) -> float:
    return float(np.mean(np.sum((points - points.mean(axis=0)) ** 2, axis=1)))


def _replicate_geography(tree, observed_tips, sigma2, error, rng):
    positions = np.zeros((tree.n_nodes, 2))
    branches = tree.branch_lengths()
    children = tree.children()
    order = [tree.root]
    points = []
    for node in order:
        order.extend(children[node])
        parent = tree.parent[node]
        if parent >= 0:
            positions[node] = positions[parent] + rng.normal(0.0, np.sqrt(sigma2 * branches[node]), 2)
        if tree.taxon[node] in observed_tips:
            points.append(positions[node] + rng.normal(0.0, error, 2) if error else positions[node])
    return _spatial_dispersion(np.array(points))


def run_posterior_predictive(bundle_dir: str, n_draws: int = 40, seed: int = 0) -> dict:
    bundle = PosteriorBundle.load(bundle_dir)
    if bundle.context.get("schema_version") != 2:
        raise ValueError("posterior predictive checks require exact study context")
    model = bundle.context["model"]
    columns = model["columns"]
    if not columns or not 1 <= n_draws <= 1000:
        raise ValueError("posterior predictive checks require nonempty columns and a bounded draw count")
    config = model["config"]
    minimum = model["minimum_present"] if config["ascertainment"] else 0
    rates = [1.0] if config["rate_alpha"] is None else _core.discrete_gamma_rates(
        config["rate_alpha"], config["rate_categories"],
    )
    mask = np.asarray(columns) >= 0
    observed = stats.all_statistics(columns)
    observed.pop("n_characters")
    coords = model["tip_coords"]
    spatial = config["use_geography"] and model["region_model"] is None and len(coords) >= 2
    if spatial:
        observed["planar_dispersion_degree2"] = _spatial_dispersion(np.array(list(coords.values())))

    rng = np.random.default_rng(seed)
    all_draws = [draw for chain in bundle.chains for draw in chain]
    if not all_draws:
        raise ValueError("empty posterior bundle")
    picks = rng.choice(len(all_draws), size=min(n_draws, len(all_draws)), replace=False)
    predicted = {name: [] for name in observed}
    for idx in picks:
        draw = all_draws[idx]
        tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"])
        replicate = _replicate(
            tree, draw["gain"] * draw["clock_rate"], draw["loss"] * draw["clock_rate"],
            rates, mask, minimum, rng, model["borrow_candidates"], draw["borrow_rates"],
            draw.get("clock_multipliers", ()),
        )
        summary = stats.all_statistics(replicate)
        if spatial:
            summary["planar_dispersion_degree2"] = _replicate_geography(
                tree, {int(t) for t in coords}, draw["geo_sigma2"],
                config.get("geography_error_degrees", 0.0), rng,
            )
        for name in predicted:
            predicted[name].append(summary[name])

    return {
        "study": bundle.dataset_name,
        "model": bundle.model_name,
        "source_commit": bundle.source_commit,
        "seed": seed,
        "draw_indices": picks.tolist(),
        "statistics": {
            name: {
                "observed": observed[name],
                "predicted_mean": float(np.mean(values)),
                "predicted_sd": float(np.std(values)),
                "bayesian_p": (
                    None if name == "missingness_rate" or (name == "singleton_fraction" and minimum == 2)
                    else float(np.mean(np.asarray(values) >= observed[name]))
                ),
                "replicates": values,
            }
            for name, values in predicted.items()
        },
        "missingness": "Conditioned on the original observed-tip masks.",
        "borrowing": "Loan annotation outcomes are outside this predictive model.",
    }
