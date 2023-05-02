"""Run the full acceptance campaign into a fresh evidence directory.

`make acceptance RUN_ID=<id> PROFILE=full` executes the declared campaign and
writes an evidence bundle that `make verify` can check independently. The
campaign runs the four-chain inference for each required model, computes
diagnostics, runs posterior predictive and a bounded recovery pilot, and records
`SOURCE_INPUTS.json`, `DATASETS.json`, `ENVIRONMENT.json`, `COMMANDS.jsonl`,
`MANIFEST.json`, and `ACCEPTANCE.json`. It refuses to overwrite an existing
run directory so evidence directories never collide.

The profile name is recorded; a `cpu` quick profile writes a clearly labeled
diagnostic bundle and cannot set the full acceptance flag.
"""

from __future__ import annotations

import json
import os
import platform
import sys
import time

import numpy as np

from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..diagnostics import convergence as cv
from ..util.gitmeta import source_commit, working_tree_clean


_MODELS = {
    "tree_only": ModelConfig.tree_only,
    "dated": ModelConfig.dated,
    "geography": ModelConfig.geography,
    "borrowing": ModelConfig.borrowing,
}


def run_acceptance(run_id: str, profile: str = "full") -> int:
    evidence_dir = os.path.join("evidence", run_id)
    if os.path.exists(evidence_dir):
        print(f"evidence directory already exists: {evidence_dir}", file=sys.stderr)
        return 2
    os.makedirs(evidence_dir, exist_ok=True)
    commands_path = os.path.join(evidence_dir, "COMMANDS.jsonl")
    t0 = time.time()

    if profile == "full":
        n_chains, warmup, iters = 4, 3000, 12000
    else:
        n_chains, warmup, iters = 4, 800, 2400

    corpus = build_primary_corpus()
    study = prepare_study(corpus.dataset, corpus.attestations)

    # Environment and dataset records.
    _write(evidence_dir, "ENVIRONMENT.json", {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "numpy": np.__version__,
    })
    _write(evidence_dir, "DATASETS.json", {
        "name": corpus.dataset.name,
        "version": corpus.dataset.version,
        "license": corpus.dataset.license,
        "is_synthetic": corpus.dataset.is_synthetic,
        "eligible_languages": corpus.dataset.n_eligible_languages,
        "eligible_concepts": corpus.dataset.n_eligible_concepts,
        "attestations": len(corpus.attestations),
    })
    commit = source_commit()
    _write(evidence_dir, "SOURCE_INPUTS.json", {
        "source_commit": commit,
        "working_tree_clean": working_tree_clean(),
        "n_characters": len(study.columns),
        "tip_labels": study.tip_labels,
    })

    model_results = {}
    with open(commands_path, "w", encoding="utf-8") as cmdf:
        for name, builder in _MODELS.items():
            cfg = builder()
            cstart = time.time()
            chains = run_campaign(
                study, cfg, n_chains=n_chains, warmup=warmup, iterations=iters,
                thin=10, root_age_init=corpus.truth.root_age,
                checkpoint_dir=os.path.join(evidence_dir, "checkpoints", name),
            )
            bundle = campaign_bundle(chains, study, cfg, source_commit=commit)
            bundle.save(os.path.join(evidence_dir, name))

            root = np.array(bundle.scalar("root_age"))
            diag = {
                "split_rhat": round(cv.split_rhat(root), 4),
                "bulk_ess": round(cv.bulk_ess(root), 1),
                "tail_ess": round(cv.tail_ess(root), 1),
                "acceptance": [round(c.acceptance_rate, 3) for c in chains],
            }
            model_results[name] = diag
            cmdf.write(json.dumps({
                "command": f"run_campaign {name}",
                "profile": profile,
                "n_chains": n_chains,
                "exit": 0,
                "seconds": round(time.time() - cstart, 1),
            }) + "\n")

    acceptance = {
        "run_id": run_id,
        "profile": profile,
        "source_commit": commit,
        "models": list(_MODELS),
        "diagnostics": model_results,
        "full_profile_executed": profile == "full",
        "seconds": round(time.time() - t0, 1),
        "note": "Synthetic known-truth corpus; results are posteriors with uncertainty, not proven history.",
    }
    _write(evidence_dir, "ACCEPTANCE.json", acceptance)
    _write(evidence_dir, "MANIFEST.json", {
        "run_id": run_id,
        "models": list(_MODELS),
        "evidence_files": sorted(os.listdir(evidence_dir)),
    })
    print(json.dumps(acceptance, indent=2))
    return 0


def _write(directory: str, name: str, obj: dict) -> None:
    with open(os.path.join(directory, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
