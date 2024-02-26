"""Evaluate declared histories across source-bound saved posterior bundles."""

from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
import re

from ..infer.artifacts import PosteriorBundle, _digest
from ..models.constraints import clade_masks
from ..models.tree import PhyloTree
from ..verify.verifier import _read_json
from .scenarios import Scenario, ScenarioRegistry, _identity


def _definitions(plan):
    if plan["schema"] != 1 or not 1 <= len(plan["scenarios"]) <= 32:
        raise ValueError("scenario plan requires one to 32 declarations")
    identifiers = set()
    for event in plan["scenarios"]:
        key = event["id"]
        if not re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", key) or key in identifiers:
            raise ValueError("scenario identifiers must be unique and stable")
        identifiers.add(key)
        if not isinstance(event["title"], str) or not event["title"].strip():
            raise ValueError("each scenario requires a title")
        for field in ("assumptions", "evidence_refs"):
            if (not isinstance(event[field], list) or not 1 <= len(event[field]) <= 16
                    or any(not isinstance(value, str) or not value.strip() for value in event[field])):
                raise ValueError(f"each scenario requires declared {field}")
        if event["kind"] == "root_age_interval":
            low, high = event["age_years"]
            if not 0 <= low < high or not math.isfinite(high) or type(event["calendar_epoch"]) is not int:
                raise ValueError("scenario age interval and epoch are invalid")
        elif event["kind"] == "clade":
            tips = event["taxa"]
            if (not isinstance(tips, list) or not 2 <= len(tips) <= 1000 or len(tips) != len(set(tips))
                    or any(not isinstance(tip, str) or not tip for tip in tips)):
                raise ValueError("clade scenarios require distinct source tip identifiers")
        else:
            raise ValueError(f"unknown scenario event {event['kind']}")
    return plan["scenarios"]


def _load_bundle(directory, source_commit):
    payload = _read_json(directory / "bundle.json", 128 * 1024 * 1024)
    bundle = PosteriorBundle(**payload)
    manifest = _read_json(directory / "MANIFEST.json", 1024 * 1024)
    if (bundle.source_commit != source_commit or manifest["source_commit"] != source_commit
            or manifest["draw_digest"] != bundle.draw_digest() or manifest["context_sha256"] != _digest(bundle.context)
            or not 1 <= bundle.n_chains == len(bundle.chains) <= 8 or any(not chain for chain in bundle.chains)
            or not 2 <= len(bundle.tip_labels) <= 1000 or len(bundle.tip_labels) != len(set(bundle.tip_labels))):
        raise ValueError("scenario bundle identity or trajectory shape differs")
    return bundle, manifest


def _evaluate(bundle, events):
    positions = {name: taxon for taxon, name in enumerate(bundle.tip_labels)}
    tip_indices = set(positions.values())
    all_tips = (1 << len(positions)) - 1
    clades = {}
    for event in events:
        if event["kind"] == "clade":
            if not set(event["taxa"]) < set(positions):
                raise ValueError(f"clade scenario does not define a proper group in {bundle.dataset_name}")
            clades[event["id"]] = sum(1 << positions[name] for name in event["taxa"])
        elif event["calendar_epoch"] != bundle.context.get("calendar_epoch"):
            raise ValueError("scenario calendar epoch differs from the bundle")

    counts = {event["id"]: [] for event in events}
    for chain in bundle.chains:
        chain_counts = dict.fromkeys(counts, 0)
        for draw in chain:
            masks = set()
            taxa, parents, heights = (draw[field] for field in ("taxon", "parent", "heights"))
            tips = [taxon for taxon in taxa if taxon >= 0]
            if (not len(parents) == len(taxa) == len(heights) <= 2 * len(positions) - 1
                    or len(tips) != len(positions) or set(tips) != tip_indices
                    or parents.count(-1) != 1 or any(parent < -1 or parent >= len(taxa) for parent in parents)):
                raise ValueError("scenario draw does not cover the declared taxa")
            age = heights[parents.index(-1)]
            if not math.isfinite(age) or age < 0 or age != draw["root_age"]:
                raise ValueError("scenario draw has an invalid root age")
            if clades:
                masks = clade_masks(PhyloTree(parents, heights, taxa))
                if all_tips not in masks:
                    raise ValueError("scenario draw does not cover the declared taxa")
            for event in events:
                key = event["id"]
                if key in clades:
                    chain_counts[key] += clades[key] in masks
                else:
                    low, high = event["age_years"]
                    chain_counts[key] += low <= age < high
        for key, count in chain_counts.items():
            counts[key].append({"matching_draws": count, "total_draws": len(chain)})
    return counts


def compare_scenarios(plan_path: Path, bundles: dict[str, Path]) -> dict:
    """Replay events using exact clade membership and half-open calendar windows."""
    plan = _read_json(plan_path, 128 * 1024)
    events = _definitions(plan)
    commits = plan["source_commits"]
    if (not 1 <= len(bundles) <= 8 or set(bundles) != set(commits)
            or any(not re.fullmatch(r"[0-9a-f]{40}", commit) for commit in commits.values())):
        raise ValueError("every comparison model requires one trusted source commit")
    registry = ScenarioRegistry()
    study = None
    for alias, directory in bundles.items():
        bundle, manifest = _load_bundle(directory, commits[alias])
        identity = (bundle.dataset_name, bundle.dataset_version, bundle.is_synthetic, bundle.tip_labels,
                    bundle.context.get("frozen_study"))
        if study is not None and study != identity:
            raise ValueError("scenario comparisons require the same frozen study and taxon order")
        study = identity
        counts = _evaluate(bundle, events)
        metadata = _identity(bundle)
        for event in events:
            parts = counts[event["id"]]
            matches = sum(row["matching_draws"] for row in parts)
            total = sum(row["total_draws"] for row in parts)
            registry.add(Scenario(
                name=event["id"], description=event["title"], model_name=bundle.model_name,
                assumptions=event["assumptions"], evidence_refs=event["evidence_refs"],
                support=matches / total,
                support_kind="clade_frequency" if event["kind"] == "clade" else "scalar_interval",
                detail={**metadata, "model_alias": alias, "context_sha256": manifest["context_sha256"],
                        "definition": event, "matching_draws": matches, "per_chain": parts},
            ))
    return {"schema": 1, "plan_sha256": _digest(plan), "scenarios": [asdict(event) for event in registry.scenarios]}
