"""Export a bounded workbench snapshot from exact saved posterior context."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from ..infer.artifacts import PosteriorBundle
from ..models.posterior import planar_geography
from ..models.tree import PhyloTree
from .sampling import run_diagnose

_MAX_DRAWS = 120
_MAX_MODELS = 8


def _retained_subset(bundle: PosteriorBundle) -> list[list[dict]]:
    if not bundle.chains or len(bundle.chains) > _MAX_DRAWS or any(not c for c in bundle.chains):
        raise ValueError("workbench requires nonempty, bounded posterior chains")
    per_chain = _MAX_DRAWS // len(bundle.chains)
    return [
        [chain[i] for i in np.linspace(0, len(chain) - 1, min(per_chain, len(chain)), dtype=int)]
        for chain in bundle.chains
    ]


def _root_locations(bundle: PosteriorBundle, chains: list[list[dict]]) -> list[dict]:
    model = bundle.context["model"]
    if not model["config"]["use_geography"] or model["region_model"] is not None:
        return []
    coords = {int(taxon): tuple(point) for taxon, point in model["tip_coords"].items()}
    if len(coords) < 2:
        return []
    rng = np.random.default_rng(20261104)
    cloud = []
    for chain in chains:
        for draw in chain:
            tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"])
            _, root = planar_geography(tree, coords, draw["geo_sigma2"])
            if root is None:
                continue
            lat, lon, variance = root
            sample = rng.normal((lat, lon), np.sqrt(variance))
            cloud.append({"lat": float(sample[0]), "lon": float(sample[1]), "age": draw["root_age"]})
    return cloud


def export_workbench(
    bundle_dirs: dict[str, str], out_path: str, source: Path | None = None,
) -> dict:
    """Export saved draws and linked data without initiating an inference chain."""
    if not 1 <= len(bundle_dirs) <= _MAX_MODELS:
        raise ValueError("workbench requires one to eight saved bundles")
    bundles = {}
    first = None
    calendar_epoch = None
    for name, directory in bundle_dirs.items():
        bundle = PosteriorBundle.load(directory)
        if bundle.context.get("schema_version") != 2:
            raise ValueError("bundle requires replayable model and evidence context")
        if first is None:
            first = bundle
            calendar_epoch = bundle.context.get(
                "calendar_epoch", 2000 if bundle.dataset_name == "iecor" else None,
            )
        if (
            bundle.tip_labels != first.tip_labels or bundle.dataset_name != first.dataset_name
            or bundle.dataset_version != first.dataset_version
            or bundle.is_synthetic != first.is_synthetic
            or bundle.context["attestations"] != first.context["attestations"]
            or bundle.context["model"]["tip_coords"] != first.context["model"]["tip_coords"]
            or bundle.context["model"]["columns"] != first.context["model"]["columns"]
        ):
            raise ValueError("workbench model bundles must share the same study observations")

        retained = _retained_subset(bundle)
        diagnostic = run_diagnose(directory)
        bundles[name] = {
            "model_name": bundle.model_name,
            "dataset_name": bundle.dataset_name,
            "dataset_version": bundle.dataset_version,
            "is_synthetic": bundle.is_synthetic,
            "n_chains": bundle.n_chains,
            "draws_per_chain": bundle.draws_per_chain,
            "tip_labels": bundle.tip_labels,
            "source_commit": bundle.source_commit,
            "notes": bundle.notes,
            "converged": diagnostic["converged"],
            "rootAgeDiagnostic": diagnostic["scalars"]["root_age"],
            "chains": [
                [{k: d[k] for k in ("iteration", "root_age", "gain", "loss", "clock_rate",
                                     "log_posterior", "newick")} for d in chain]
                for chain in retained
            ],
            "rootLocations": _root_locations(bundle, retained),
            "calendarEpoch": calendar_epoch,
        }

    names = first.context.get("language_names", {})
    if source is not None:
        from ..io.freeze import verify_inputs

        inputs = first.context["frozen_study"]["source"]["files"]
        filename = "cldf/languages.csv"
        verify_inputs(source, {filename: inputs[filename]})
        with (source / filename).open(newline="", encoding="utf-8") as raw:
            names = {r["ID"]: r["Name"] for r in csv.DictReader(raw)}
    names = {tip: names.get(tip, tip) for tip in first.tip_labels}
    for bundle in bundles.values():
        bundle["tip_names"] = names

    coordinates = first.context["model"]["tip_coords"]
    modern_points = [
        {"lat": lat, "lon": lon, "kind": "modern-reference", "label": first.tip_labels[int(t)]}
        for t, (lat, lon) in coordinates.items()
    ]
    attestations = []
    for a in first.context["attestations"]:
        if a["language_id"] not in first.tip_labels:
            continue
        # EDH age values in frozen bundles use the source cohort's 1950 epoch.
        # IE-CoR dates use 2000 CE. Convert only for the linked display.
        source_epoch = 1950 if a["id"].startswith("edh:") else None
        shift = calendar_epoch - source_epoch if calendar_epoch and source_epoch else 0
        attestations.append({
            "id": a["id"], "lat": a["latitude"], "lon": a["longitude"], "kind": "attestation",
            "label": a["language_id"], "ageLow": a["age_low"] + shift, "ageHigh": a["age_high"] + shift,
            "sourceAgeLow": a["age_low"], "sourceAgeHigh": a["age_high"], "sourceEpoch": source_epoch,
            "uncertaintyKm": a["location_uncertainty_km"], "region": a["region"],
            "sourceRef": a["source_ref"], "page": a["page"], "bibliography": a["bibliography"],
            "locationKind": a["location_kind"], "locationNote": a["location_uncertainty_note"],
            "languageNote": a["language_assignment_note"], "conflicts": a["conflicting_with"],
        })
    matrix = {
        "languageLabels": [names[tip] for tip in first.tip_labels],
        "columns": first.context["model"]["columns"][:400],
    }
    payload = {
        "bundles": bundles, "modernPoints": modern_points,
        "attestations": attestations, "matrix": matrix,
    }
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    return {
        "out": str(out), "models": list(bundles), "modern_points": len(modern_points),
        "attestations": len(attestations), "matrix_columns": len(matrix["columns"]),
    }
