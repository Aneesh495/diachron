"""Run and record the bounded diagnostic campaign.

Full acceptance fails closed until the observed-data prediction, sensitivity,
and release gates can be executed and independently verified.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from .sampling import run_diagnose
from ..util.gitmeta import source_commit, working_tree_clean


_MODELS = {
    "tree_only": ModelConfig.tree_only,
    "dated": ModelConfig.dated,
    "geography": ModelConfig.geography,
    "borrowing": ModelConfig.borrowing,
}


def run_acceptance(run_id: str, profile: str = "full") -> int:
    if profile not in {"cpu", "full"} or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", run_id):
        print("invalid profile or run ID", file=sys.stderr)
        return 2
    evidence_dir = os.path.join("evidence", run_id)
    try:
        os.makedirs(evidence_dir, exist_ok=False)
    except FileExistsError:
        print(f"evidence directory already exists: {evidence_dir}", file=sys.stderr)
        return 2
    commit = source_commit()
    if profile == "full":
        result = {
            "run_id": run_id, "profile": profile, "source_commit": commit,
            "status": "blocked", "full_profile_executed": False,
            "blocking_gates": ["G06", "G07", "G08", "G10", "G14", "G15"],
            "reason": "Full acceptance needs the complete observed-data campaign and verifier.",
        }
        _write(evidence_dir, "ACCEPTANCE.json", result)
        print(json.dumps(result, indent=2))
        return 1

    commands_path = os.path.join(evidence_dir, "COMMANDS.jsonl")
    t0 = time.time()
    n_chains, warmup, iters = 4, 400, 1200
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
    tracked = subprocess.check_output(["git", "ls-files", "-z"]).decode().split("\0")
    inputs = {
        name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
        for name in tracked if name and Path(name).is_file()
    }
    _write(evidence_dir, "SOURCE_INPUTS.json", {
        "source_commit": commit,
        "working_tree_clean": working_tree_clean(),
        "n_characters": len(study.columns),
        "tip_labels": study.tip_labels,
        "tracked_file_sha256": inputs,
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

            diag = run_diagnose(os.path.join(evidence_dir, name))
            diag["acceptance"] = [round(c.acceptance_rate, 3) for c in chains]
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
        "status": "diagnostic_completed",
        "full_profile_executed": False,
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
        json.dump(obj, f, indent=2, allow_nan=False)
