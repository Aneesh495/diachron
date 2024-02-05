"""Known-truth recovery study (gate G05).

We generate many synthetic datasets spanning tree shape, missingness, borrowing,
and dating regimes, fit a short chain to each, and measure two things: whether the
inferred tree recovers the true clades (clade recovery rate) and whether the true
root age falls inside the posterior credible interval (credible-region coverage).

Coverage is the honesty check: if the 90% credible interval contains the truth
about 90% of the time across datasets, the uncertainty is calibrated. Systematic
over- or under-coverage is reported rather than hidden, because a posterior that
is narrower than it should be is a failure even when its point estimate is close.

Each dataset is labeled synthetic. The study runs a bounded number of short
chains; the count and per-dataset budget are arguments so a pilot and the full
200-dataset run use the same code path.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

from ..io.simulate import simulate_dataset, simulate_attestations
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..diagnostics import topology as tp
from ..models.tree import tree_from_truth
from ..util.gitmeta import source_commit


def _true_clades(truth) -> set:
    tree = tree_from_truth(truth)
    nw = tree.newick()
    return {clade for clade in tp._clade_signature(nw) if len(clade) < len(truth.tip_language_ids)}


def _regime(i: int) -> dict:
    """Vary the generative regime across datasets to span the required space."""
    rng = np.random.default_rng(10_000 + i)
    return {
        "n_languages": int(rng.integers(8, 20)),
        "n_concepts": int(rng.integers(40, 90)),
        "root_age": float(rng.uniform(3000, 8000)),
        "gain": float(rng.uniform(0.0002, 0.0006)),
        "loss": float(rng.uniform(0.0006, 0.0015)),
        "borrow_prob": float(rng.choice([0.0, 0.15, 0.3])),
        "missing_frac": float(rng.choice([0.0, 0.08, 0.2])),
    }


def run_recovery_study(n_datasets: int = 200, out_dir: str = "runs/recovery",
                       warmup: int = 400, iterations: int = 1200) -> dict:
    if not 1 <= n_datasets <= 5000:
        raise ValueError("recovery study requires between 1 and 5000 datasets")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    per_dataset = []
    commit = source_commit()
    with (out / "rows.jsonl").open("w", encoding="utf-8") as raw:
        for i in range(n_datasets):
            reg = _regime(i)
            data_seed = 20000 + i
            chain_seed = 30000 + i
            row = {"i": i, "regime": reg, "data_seed": data_seed, "chain_seed": chain_seed}
            try:
                dataset, truth = simulate_dataset(seed=data_seed, name=f"rec{i}", **reg)
                attestations = simulate_attestations(
                    dataset, truth, n_records=(12, 20, 32)[i % 3], seed=data_seed,
                )
                study = prepare_study(dataset, attestations)
                calibration = study.calibrations[0]
                initial = (calibration.age_low + calibration.age_high) / 2.0
                config = ModelConfig.dated()
                chains = run_campaign(
                    study, config, n_chains=2, warmup=warmup, iterations=iterations, thin=10,
                    root_age_init=initial, base_seed=chain_seed, initialization="random",
                )
                bundle = campaign_bundle(chains, study, config, source_commit=commit)
                newicks = [nw for c in range(bundle.n_chains) for nw in bundle.chain_newicks(c)]
                post_clades = tp.clade_posteriors(newicks)
                truth_clades = _true_clades(truth)
                masses = [
                    post_clades.get(clade, 0.0) for clade in sorted(
                        truth_clades, key=lambda taxa: tuple(sorted(taxa)),
                    )
                ]
                root_samples = [d for ch in bundle.scalar("root_age") for d in ch]
                lo, hi = np.quantile(root_samples, [0.05, 0.95])
                row.update({
                    "status": "complete", "n_attestations": len(attestations),
                    "root_true": truth.root_age, "root_ci": [float(lo), float(hi)],
                    "root_in_90ci": bool(lo <= truth.root_age <= hi),
                    "clade_recovery": float(np.mean(np.array(masses) > 0.5)),
                    "root_samples": root_samples, "true_clade_masses": masses,
                })
            except (ValueError, ArithmeticError) as exc:
                row.update({"status": "failed", "error": str(exc)})
            raw.write(json.dumps(row, allow_nan=False) + "\n")
            raw.flush()
            per_dataset.append(row)

    completed = [row for row in per_dataset if row["status"] == "complete"]
    result = {
        "gate": "G05", "source_commit": commit,
        "profile": "known-truth recovery, protocol 2",
        "n_datasets": n_datasets, "completed": len(completed),
        "failed": n_datasets - len(completed),
        "iterations": iterations, "warmup": warmup, "chains": 2, "thin": 10,
        "mean_clade_recovery": float(np.mean([r["clade_recovery"] for r in completed]))
        if completed else None,
        "root_age_90ci_coverage": float(np.mean([r["root_in_90ci"] for r in completed]))
        if completed else None,
        "coverage_target": 0.90,
    }
    compressed = out / "rows.jsonl.gz"
    with compressed.open("wb") as output, (out / "rows.jsonl").open("rb") as source:
        with gzip.GzipFile(filename="", fileobj=output, mode="wb", mtime=0) as archive:
            shutil.copyfileobj(source, archive)
    with compressed.open("rb") as source:
        result["raw_rows_sha256"] = hashlib.file_digest(source, "sha256").hexdigest()
    result["covered"] = sum(row["root_in_90ci"] for row in completed)
    groups = defaultdict(list)
    for row in completed:
        groups[str(row["regime"]["borrow_prob"])].append(row)
    result["by_borrow_probability"] = {
        key: {"datasets": len(rows), "coverage": float(np.mean([row["root_in_90ci"] for row in rows])),
              "mean_clade_recovery": float(np.mean([row["clade_recovery"] for row in rows]))}
        for key, rows in sorted(groups.items())
    }
    (out / "recovery-summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (out / "recovery.json").write_text(
        json.dumps({"summary": result, "per_dataset": per_dataset}, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    return result
