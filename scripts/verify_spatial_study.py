"""Recompute the study scores and replay selected emissions with SciPy expm."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import fmean

import numpy as np
from scipy.linalg import expm

from diachron.geo.study import _AXES, _case, _generator, fibonacci_sphere


def close(actual, expected):
    if not np.allclose(actual, expected, rtol=0, atol=1e-9):
        raise ValueError(f"spatial replay differs: {actual!r} versus {expected!r}")


def vector(values):
    if len(values) != 6 or any(not 0 <= p <= 1 for p in values):
        raise ValueError("invalid spatial predictive probabilities")
    close(sum(values), 1.)
    return values


def replay(row, truth, models, scales):
    case = _case(row["seed"], truth, models["regions"][0])
    tree = case.tree
    if (row["root_region"], row["tip_regions"], row["withheld"], row["parent"], row["heights"]) != (
        case.root_region, case.tip_regions, case.withheld, tree.parent, tree.heights,
    ):
        raise ValueError("known spatial truth differs from its frozen seed")
    children = tree.children()
    nodes = [tree.root]
    for node in nodes:
        nodes.extend(children[node])
    for name, (generator, mapping) in models.items():
        transitions = [expm(generator * (scales[name] * time)) for time in tree.branch_lengths()]
        observed = {tree.taxon.index(taxon): region for taxon, region in enumerate(case.tip_regions)
                    if taxon not in case.withheld}

        def likelihood(evidence):
            messages = {}
            for node in reversed(nodes):
                message = np.ones(len(mapping)) if node not in evidence else (mapping == evidence[node]).astype(float)
                for child in children[node]:
                    message *= transitions[child] @ messages[child]
                messages[node] = message
            return float(messages[tree.root].mean())

        denominator = likelihood(observed)
        for node, expected in zip([tree.root, *case.withheld], [row["predictions"][name]["root"],
                                                               *row["predictions"][name]["heldout"]]):
            close([likelihood(observed | {node: region}) / denominator for region in range(6)], expected)


def verify(summary_path, cases_path, expected_commit):
    if summary_path.stat().st_size > 65536 or cases_path.stat().st_size > 2_000_000:
        raise ValueError("oversized spatial study artifacts")
    summary = json.loads(summary_path.read_text())
    if (summary["source_commit"] != expected_commit or summary["schema"] != 1 or summary["is_synthetic"] is not True
            or [summary[name] for name in ("train_per_truth", "test_per_truth", "regions", "mesh_states", "test_tips_per_tree")]
            != [40, 120, 6, 96, 3]):
        raise ValueError("spatial study is outside the frozen protocol")
    content = cases_path.read_bytes()
    rows = [json.loads(line) for line in content.splitlines()]
    if len(rows) != 240:
        raise ValueError("spatial study requires all 240 test cases")
    mesh = fibonacci_sphere(96)
    models = {"regions": (_generator(_AXES, local=False), np.arange(6)),
              "sphere_mesh": (_generator(mesh, local=True), np.argmax(mesh @ _AXES.T, axis=1))}
    for index, truth in enumerate(("sphere", "regions")):
        result = summary["results"][truth]
        scores = {name: [] for name in models}
        root_sizes = {name: [] for name in models}
        coverage = {name: [] for name in models}
        for offset, row in enumerate(rows[120 * index:120 * (index + 1)]):
            if row["truth"] != truth or row["seed"] != 64200 + 1000 * index + offset or len(set(row["withheld"])) != 3:
                raise ValueError("wrong spatial test partition")
            for name in models:
                prediction = row["predictions"][name]
                if len(prediction["heldout"]) != 3:
                    raise ValueError("missing withheld tip prediction")
                scores[name].append(-fmean(math.log(vector(p)[row["tip_regions"][tip]])
                                           for p, tip in zip(prediction["heldout"], row["withheld"])))
                root = vector(prediction["root"])
                chosen = []
                for region in sorted(range(6), key=lambda region: -root[region]):
                    chosen.append(region)
                    if sum(root[state] for state in chosen) >= .9:
                        break
                covered = row["root_region"] in chosen
                if prediction["root_90set"] != {"size": len(chosen), "covered": covered}:
                    raise ValueError("incorrect root credible set")
                root_sizes[name].append(len(chosen))
                coverage[name].append(covered)
            if offset in {0, 39, 79, 119}:
                replay(row, truth, models, result["selected_rate_scales"])
        for name in models:
            score = result["scores"][name]
            close([fmean(scores[name]), fmean(coverage[name]), fmean(root_sizes[name])],
                  [score["heldout_nll"], score["root_90set_coverage"], score["root_90set_size"]])
            training = result["training_log_likelihood"][name]
            close(float(max(training, key=training.get)), result["selected_rate_scales"][name])
        difference = np.subtract(scores["sphere_mesh"], scores["regions"])
        rng = np.random.default_rng(65327)
        boot = difference[rng.integers(0, 120, (1000, 120))].mean(axis=1)
        paired = result["scores"]["paired_nll_mesh_minus_regions"]
        close(paired["mean"], fmean(difference))
        close(paired["bootstrap_95ci"], np.quantile(boot, [.025, .975]))
    return {"passed": True, "test_cases": len(rows), "independent_transition_probes": 8,
            "cases_sha256": hashlib.sha256(content).hexdigest(), "source_commit": expected_commit}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--expected-commit", required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.summary, args.cases, args.expected_commit), indent=2))
