"""Export a workbench.json the TypeScript UI can load (bounded serving, D35).

This reads already-saved posterior bundles and the dataset context and writes a
single JSON file with everything the four linked views need: the bundles (one per
model), the modern reference coordinates, the dated attestations, and the
character matrix. It performs no inference; moving the timeline or switching
models in the UI only reads this file, so a user interaction can never trigger an
uncontrolled chain.

Geographic honesty is enforced at the boundary: modern coordinates are tagged
`modern-reference`, attestations `attestation`, and the UI renders inferred
latent locations as clouds from the posterior draws. The exporter never writes a
single crisp homeland coordinate.
"""

from __future__ import annotations

import json
import os

from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study
from ..infer.artifacts import PosteriorBundle


def export_workbench(bundle_dirs: dict[str, str], out_path: str) -> dict:
    """Write workbench.json from a mapping of model name -> bundle directory."""
    corpus = build_primary_corpus()
    study = prepare_study(corpus.dataset, corpus.attestations)

    bundles = {}
    for name, d in bundle_dirs.items():
        if not os.path.isdir(d):
            continue
        b = PosteriorBundle.load(d)
        bundles[name] = {
            "model_name": b.model_name,
            "dataset_name": b.dataset_name,
            "dataset_version": b.dataset_version,
            "is_synthetic": b.is_synthetic,
            "n_chains": b.n_chains,
            "draws_per_chain": b.draws_per_chain,
            "tip_labels": b.tip_labels,
            "chains": b.chains,
            "source_commit": b.source_commit,
            "notes": b.notes,
        }

    # Modern reference coordinates in the tip-label order the bundle uses.
    coord_by_lang = {
        lang.id: (lang.latitude, lang.longitude) for lang in corpus.dataset.languages
    }
    modern_points = []
    for lid in study.tip_labels:
        lat, lon = coord_by_lang.get(lid, (None, None))
        modern_points.append(
            {
                "lat": lat if lat is not None else 0.0,
                "lon": lon if lon is not None else 0.0,
                "kind": "modern-reference",
                "label": lid,
            }
        )

    attestations = [
        {
            "lat": a.latitude if a.latitude is not None else 0.0,
            "lon": a.longitude if a.longitude is not None else 0.0,
            "kind": "attestation",
            "label": a.language_id,
            "uncertaintyKm": a.location_uncertainty_km,
            "ageLow": a.age_low,
            "ageHigh": a.age_high,
        }
        for a in corpus.attestations
        if a.language_id in study.tip_labels
    ]

    matrix = {
        "languageLabels": study.tip_labels,
        "columns": study.columns[:400],  # cap columns for responsive rendering
    }

    payload = {
        "bundles": bundles,
        "modernPoints": modern_points,
        "attestations": attestations,
        "matrix": matrix,
    }
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f)
    return {
        "out": out_path,
        "models": list(bundles),
        "modern_points": len(modern_points),
        "attestations": len(attestations),
        "matrix_columns": len(matrix["columns"]),
    }
