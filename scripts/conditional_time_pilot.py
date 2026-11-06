"""Isolate within-topology calendar mixing on two saved observed trees."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from functools import partial
import json
import multiprocessing
from pathlib import Path
import sys
import time

from diachron.cli.sampling import run_diagnose
from diachron.infer.artifacts import bundle_from_chains
from diachron.infer.chain import Chain, ChainConfig, ChainState
from diachron.models.constraints import clade_masks
from diachron.models.tree import PhyloTree
from diachron.registry.compare import _load_bundle
from diachron.util.gitmeta import source_commit, working_tree_clean
from diachron.verify.verifier import _model, verify_evidence


def fit(index, *, model, anchor, labels, destination):
    tree = PhyloTree(list(anchor["parent"]), list(anchor["heights"]), list(anchor["taxon"]), dict(enumerate(labels)))
    children = tree.children()
    base = max(tree.heights[child] for child in children[tree.root])
    tree.heights[tree.root] = base + (tree.root_age() - base) * (.75, .9, 1.1, 1.3)[index]
    state = ChainState(anchor["gain"] * (.7, .9, 1.1, 1.4)[index], anchor["loss"],
                       anchor["clock_rate"] * (1.4, 1.1, .9, .7)[index])
    config = ChainConfig(seed=84000 + 7919 * index, warmup=1000, iterations=3000, thin=10,
                         initialization="saved-conditioned-tree")
    chain = Chain(model, tree, state, config)
    path = destination / f"chain-{index}.json"
    chain.run(checkpoint_path=path, checkpoint_every=1000)
    return str(path)


def run(directory: Path, destination: Path, expected: str):
    if not working_tree_clean():
        raise ValueError("conditional-time pilot requires a clean source checkout")
    bundle, manifest = _load_bundle(directory, expected)
    if bundle.model_name != "dated" or bundle.is_synthetic or bundle.n_chains != 4:
        raise ValueError("conditional-time pilot requires the observed dated reference")
    destination.mkdir(parents=True)
    commit = source_commit()
    report = {"source_commit": commit, "anchor_source_commit": expected,
              "anchor_draw_digest": manifest["draw_digest"], "topologies": {}}
    for chain_id in (0, 3):
        started = time.monotonic()
        label = f"anchor-chain-{chain_id}"
        out = destination / label
        out.mkdir()
        anchor = bundle.chains[chain_id][0]
        tree = PhyloTree(anchor["parent"], anchor["heights"], anchor["taxon"])
        full_clade = (1 << len(bundle.tip_labels)) - 1
        constraints = [[tip for tip in range(len(bundle.tip_labels)) if mask & (1 << tip)]
                       for mask in sorted(clade_masks(tree) - {full_clade})]
        raw_model = {**bundle.context["model"], "required_clades": constraints}
        model = _model(raw_model)
        worker = partial(fit, model=model, anchor=anchor, labels=bundle.tip_labels, destination=out)
        with ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("spawn")) as pool:
            paths = list(pool.map(worker, range(4)))
        chains = [Chain.resume(model, path) for path in paths]
        fitted = bundle_from_chains(chains, bundle.model_name, bundle.dataset_name, bundle.dataset_version,
                                    False, bundle.tip_labels, commit)
        fitted.context = {**bundle.context, "model": asdict(model),
                          "chain_configurations": [asdict(chain.config) for chain in chains],
                          "move_statistics": [{"proposed": chain.move_proposals, "accepted": chain.move_accepts}
                                              for chain in chains],
                          "conditional_time": {"anchor_chain": chain_id, "anchor_iteration": anchor["iteration"],
                                               "anchor_source_commit": expected}}
        fitted.save(out / "dated")
        diagnosis = run_diagnose(out / "dated")
        if any(row["distinct"] != 1 for row in diagnosis["topology_per_chain"]):
            raise ValueError("conditional-time chain left the frozen topology")
        verification = verify_evidence(str(out / "dated"), expected_commit=commit)
        if not verification["ok"]:
            raise ValueError("conditional-time trajectories did not replay")
        report["topologies"][label] = {"anchor_iteration": anchor["iteration"], "anchor_root_age": anchor["root_age"],
            "constrained_clades": len(constraints), "seconds": time.monotonic() - started,
            "draw_digest": fitted.draw_digest(), "scalars": diagnosis["scalars"], "verified": True}
    (destination / "conditional-time.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
