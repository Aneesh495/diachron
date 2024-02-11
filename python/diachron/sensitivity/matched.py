"""Known-truth recovery under the same retained-pattern model used for fitting."""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from functools import partial
import fcntl
import gzip
import hashlib
import json
import math
import multiprocessing
from pathlib import Path

import numpy as np
from scipy.stats import truncnorm

from ..diagnostics import convergence, topology
from ..infer.campaign import StudyData, run_campaign
from ..io.types import Concept, Dataset, Language
from ..models.conditional import sample_selected_patterns
from ..models.posterior import ModelConfig
from ..models.tip_dates import TipDate
from ..models.tree import PhyloTree
from ..util.gitmeta import source_commit


def serial_coalescent(labels, tip_ages, population, rng):
    """Generate waiting intervals conditional on the serial sampling times."""
    n_tips = len(labels)
    if n_tips < 2 or len(tip_ages) != n_tips or not math.isfinite(population) or population <= 0:
        raise ValueError("invalid coalescent dimensions or population scale")
    if any(not math.isfinite(age) or age < 0 for age in tip_ages):
        raise ValueError("tip ages must be finite and nonnegative")
    schedule = sorted(range(n_tips), key=lambda index: tip_ages[index])
    parents = [-1] * (2 * n_tips - 1)
    ages = list(tip_ages) + [0.0] * (n_tips - 1)
    taxa = list(range(n_tips)) + [-1] * (n_tips - 1)
    active = []
    current = 0.0
    sampling_index = 0
    next_node = n_tips
    while sampling_index < n_tips or len(active) > 1:
        next_sample = tip_ages[schedule[sampling_index]] if sampling_index < n_tips else math.inf
        count = len(active)
        waiting = rng.exponential(population / math.comb(count, 2)) if count >= 2 else math.inf
        event = current + waiting
        if next_sample <= event:
            current = next_sample
            while sampling_index < n_tips and tip_ages[schedule[sampling_index]] == current:
                active.append(schedule[sampling_index])
                sampling_index += 1
        else:
            current = event
            left = active.pop(int(rng.integers(len(active))))
            right = active.pop(int(rng.integers(len(active))))
            parents[left] = parents[right] = next_node
            ages[next_node] = current
            active.append(next_node)
            next_node += 1
    return PhyloTree(parents, ages, taxa, dict(enumerate(labels)))


def generate_matched_case(index: int):
    design = np.random.default_rng(51000 + index)
    rng = np.random.default_rng(52000 + index)
    count = int(design.integers(8, 20))
    characters = int(design.integers(40, 90))
    n_dated = (0, 2, 4)[index % 3]
    missing = (.0, .08, .2)[(index // 3) % 3]
    copying = bool((index // 9) % 2)
    config = ModelConfig.borrowing() if copying else ModelConfig.dated()
    labels = [f"tip-{taxon:02d}" for taxon in range(count)]
    dated_taxa = sorted(design.choice(count, size=n_dated, replace=False).tolist())
    dates = {taxon: TipDate("normal", 700.0 + 400 * offset, 80.0 + 40 * (index % 3))
             for offset, taxon in enumerate(dated_taxa)}
    ages = [0.0] * count
    for taxon, date in dates.items():
        ages[taxon] = float(truncnorm.rvs(-date.mean / date.sd, np.inf, loc=date.mean, scale=date.sd, random_state=rng))
    tree = serial_coalescent(labels, ages, config.pop_size, rng)
    gain, loss = (float(value) for value in rng.lognormal(math.log(config.rate_prior_median), config.rate_prior_log_sd, 2))
    clock = float(rng.lognormal(0, config.clock_prior_log_sd))
    candidate_pairs = [(0, count - 1)] if copying else []
    rates = rng.exponential(1 / config.borrowing_prior_rate, 2).tolist() if copying else []
    contacts = [(donor, recipient, rate) for (donor, recipient), rate in
                zip(((0, count - 1), (count - 1, 0)), rates)]
    masks = rng.random((characters, count)) >= missing
    for mask in masks:
        if mask.sum() < 3:
            available = np.flatnonzero(~mask)
            mask[rng.choice(available, size=3 - int(mask.sum()), replace=False)] = True
    columns = sample_selected_patterns(tree, gain * clock, loss * clock, masks, 2, rng, contacts=contacts)
    dataset = Dataset(
        name=f"matched-{index}", version="3", license="CC0-1.0", doi=None,
        concepts=[Concept(str(character), f"Retained character {character}") for character in range(characters)],
        languages=[Language(label, label, is_historical=taxon in dates) for taxon, label in enumerate(labels)],
        lexemes=[], judgments=[], is_synthetic=True, source_url="synthetic:matched-binary-retention:3",
    )
    study = StudyData(labels, columns, [], {}, dataset, borrow_candidates=candidate_pairs, minimum_present=2,
                      tip_dates=dates, character_ids=[(str(character), "present") for character in range(characters)])
    regime = {"n_languages": count, "n_concepts": characters, "n_dated_tips": n_dated,
              "missing_frac": missing, "borrow_regime": "candidate_copy" if copying else "vertical",
              "root_age": tree.root_age()}
    truth = {"parent": tree.parent, "heights": tree.heights, "taxon": tree.taxon, "gain": gain,
             "loss": loss, "clock_rate": clock, "borrow_rates": rates, "tip_dates": {str(key): asdict(value)
             for key, value in dates.items()}, "columns": columns}
    return study, config, tree, regime, truth


def _case(index: int, warmup: int, iterations: int) -> dict:
    row = {"i": index, "data_seed": 52000 + index, "chain_seed": 53000 + index}
    try:
        study, config, tree, regime, truth = generate_matched_case(index)
        row.update({"regime": regime, "truth": truth})
        chains = run_campaign(study, config, n_chains=4, warmup=warmup, iterations=iterations,
                              thin=10, root_age_init=8000.0, base_seed=row["chain_seed"])
        raw_root = np.asarray([chain.scalar_trace("root_age") for chain in chains])
        roots = raw_root.flatten()
        low, high = np.quantile(roots, [.05, .95])
        true_clades = topology._clade_signature(tree.newick())
        frequencies = topology.clade_posteriors([draw.newick for chain in chains for draw in chain.draws])
        masses = [frequencies.get(group, 0.0) for group in sorted(true_clades, key=lambda taxa: tuple(sorted(taxa)))
                  if len(group) < len(study.tip_labels)]
        rhat = convergence.rank_normalized_rhat(raw_root)
        row.update({"status": "complete", "root_true": tree.root_age(), "root_ci": [float(low), float(high)],
                    "root_in_90ci": bool(low <= tree.root_age() <= high),
                    "root_rank": int(np.sum(roots < tree.root_age())), "root_samples": roots.tolist(),
                    "true_clade_masses": masses, "clade_recovery": sum(value > .5 for value in masses) / len(masses),
                    "chain_seeds": [chain.config.seed for chain in chains],
                    "root_rank_rhat": float(rhat) if np.isfinite(rhat) else None,
                    "root_bulk_ess": convergence.bulk_ess(raw_root)})
    except (ArithmeticError, ValueError) as exc:
        row.update({"status": "failed", "error": str(exc)})
    return row


def _write_json(path: Path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, allow_nan=False, indent=2) + "\n")
    temporary.replace(path)


def run_matched_recovery(n_datasets: int, out: Path, jobs: int = 4, resume: bool = False,
                         *, warmup: int = 800, iterations: int = 3200):
    if not 1 <= n_datasets <= 200 or not 1 <= jobs <= 4 or warmup < 0 or iterations < 10:
        raise ValueError("invalid matched recovery budget")
    plan = {"schema": 1, "protocol": 3, "source_commit": source_commit(), "datasets": n_datasets,
            "warmup": warmup, "iterations": iterations, "chains": 4, "thin": 10}
    unit_dir = out / "units"
    if out.exists():
        if not resume or json.loads((out / "MATCHED_PLAN.json").read_text()) != plan:
            raise ValueError("matched recovery output exists with a different source or plan")
    else:
        unit_dir.mkdir(parents=True)
        _write_json(out / "MATCHED_PLAN.json", plan)
    with (out / "MATCHED.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise FileExistsError("matched recovery is already active") from exc
        return _run_units(out, unit_dir, plan, jobs)


def _run_units(out, unit_dir, plan, jobs):
    n_datasets, warmup, iterations = plan["datasets"], plan["warmup"], plan["iterations"]
    rows = {}
    for path in sorted(unit_dir.glob("case-*.json")):
        if path.stat().st_size > 256 * 1024:
            raise ValueError("oversized saved recovery unit")
        row = json.loads(path.read_text())
        if row["source_commit"] != plan["source_commit"]:
            raise ValueError("saved recovery unit has different source")
        index = row["i"]
        if not 0 <= index < n_datasets or index in rows or path.name != f"case-{index:03d}.json":
            raise ValueError("duplicate or invalid saved recovery unit")
        rows[index] = row
    remaining = [index for index in range(n_datasets) if index not in rows]
    worker = partial(_case, warmup=warmup, iterations=iterations)

    def save(results):
        for row in results:
            row["source_commit"] = plan["source_commit"]
            _write_json(unit_dir / f"case-{row['i']:03d}.json", row)
            rows[row["i"]] = row

    if jobs == 1:
        save(map(worker, remaining))
    elif remaining:
        with ProcessPoolExecutor(max_workers=jobs, mp_context=multiprocessing.get_context("spawn")) as pool:
            save(pool.map(worker, remaining))
    ordered = [rows[index] for index in range(n_datasets)]
    compressed = out / "rows.jsonl.gz"
    with compressed.open("wb") as stream, gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as target:
        for row in ordered:
            target.write((json.dumps(row, allow_nan=False) + "\n").encode())
    completed = [row for row in ordered if row["status"] == "complete"]
    covered = sum(row["root_in_90ci"] for row in completed)
    groups = defaultdict(list)
    for row in completed:
        groups[row["regime"]["borrow_regime"]].append(row)
    histogram = [0] * 10
    for row in completed:
        bin_index = min(9, 10 * row["root_rank"] // (len(row["root_samples"]) + 1))
        histogram[bin_index] += 1
    with compressed.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    result = {"gate": "G05", "profile": "matched binary retention recovery", "protocol": 3,
              "source_commit": plan["source_commit"], "n_datasets": n_datasets,
              "completed": len(completed), "failed": n_datasets - len(completed),
              "chains": 4, "warmup": warmup, "iterations": iterations, "thin": 10,
              "covered": covered, "root_age_90ci_coverage": covered / len(completed) if completed else None,
              "coverage_target": .9, "mean_clade_recovery": float(np.mean([r["clade_recovery"] for r in completed]))
              if completed else None, "raw_rows_sha256": digest, "root_rank_histogram": histogram,
              "by_borrow_regime": {key: {"datasets": len(values),
                  "coverage": sum(row["root_in_90ci"] for row in values) / len(values),
                  "mean_clade_recovery": float(np.mean([row["clade_recovery"] for row in values]))}
                  for key, values in sorted(groups.items())}}
    _write_json(out / "recovery-summary.json", result)
    return result
