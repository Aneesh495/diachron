"""Reference corpora for the primary and second studies, plus manifests.

Diachron's experiments need lexical data at the scale the brief sets (at least
100 varieties and 100 concepts for the primary study, at least 30 for the
second). The restricted real collections (IE-CoR / Lexibank) cannot be
redistributed, so the repository ships two things instead:

1. A CLDF reader (`io.cldf`) that ingests those collections from a directory the
   user has acquired under the collection's own license, and an acquisition
   recipe in the data manifest.

2. A reproducible reference corpus built here, clearly labeled synthetic, at the
   required scale with known ground truth. Known truth is what the recovery and
   calibration gates (G05, and the coverage checks) actually need: a corpus whose
   history is known lets Diachron measure whether inference recovers it and
   whether the credible regions are honest. The synthetic flag is propagated into
   every artifact and view so it is never presented as observed history.

The corpus builder draws an Indo-European-shaped tree (a root around 6000 years
before present, major clades, and within-clade diversification) so the study has
realistic structure, then simulates cognate evolution, a documented loan set, and
dated attestations consistent with that truth.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os

from .simulate import simulate_dataset, simulate_attestations, TrueHistory
from .types import Dataset, Attestation


@dataclass
class Corpus:
    dataset: Dataset
    truth: TrueHistory
    attestations: list[Attestation]


def build_primary_corpus(seed: int = 20240601) -> Corpus:
    """Primary study: 112 varieties, 110 concepts, with loans and attestations.

    The counts exceed the 100/100 floor and leave room for held-out partitions
    without dropping below the eligible-dimension requirement.
    """
    dataset, truth = simulate_dataset(
        n_languages=112,
        n_concepts=110,
        root_age=6000.0,
        gain=0.00035,
        loss=0.00085,
        clock_rate=1.0,
        borrow_prob=0.25,     # roughly a quarter of concepts see a borrowing event
        missing_frac=0.08,    # realistic incomplete coverage
        seed=seed,
        name="diachron_primary",
    )
    attestations = simulate_attestations(dataset, truth, n_records=140, seed=seed)
    return Corpus(dataset=dataset, truth=truth, attestations=attestations)


def build_second_corpus(seed: int = 20240602) -> Corpus:
    """Second study: 36 varieties, 90 concepts, a different calibration regime.

    The second family uses its own root age and clock, so the primary study's
    calibration assumptions are not forced onto it, matching the brief's warning
    against reusing one family's calibration.
    """
    dataset, truth = simulate_dataset(
        n_languages=36,
        n_concepts=90,
        root_age=3500.0,      # a younger, shallower family
        gain=0.0005,
        loss=0.0012,
        clock_rate=1.0,
        borrow_prob=0.3,
        missing_frac=0.12,
        seed=seed,
        name="diachron_second",
    )
    attestations = simulate_attestations(dataset, truth, n_records=60, seed=seed)
    return Corpus(dataset=dataset, truth=truth, attestations=attestations)


def write_manifest(corpus: Corpus, path: str) -> dict:
    """Write a dataset manifest recording dimensions, licenses, and provenance."""
    ds = corpus.dataset
    manifest = {
        "name": ds.name,
        "version": ds.version,
        "license": ds.license,
        "doi": ds.doi,
        "is_synthetic": ds.is_synthetic,
        "notes": ds.notes,
        "dimensions": {
            "languages_total": len(ds.languages),
            "concepts_total": len(ds.concepts),
            "eligible_languages": ds.n_eligible_languages,
            "eligible_concepts": ds.n_eligible_concepts,
            "lexemes": len(ds.lexemes),
            "judgments": len(ds.judgments),
            "observed_cells": ds.observed_cell_count(),
            "attestations": len(corpus.attestations),
        },
        "true_history": {
            "root_age": corpus.truth.root_age,
            "gain": corpus.truth.gain,
            "loss": corpus.truth.loss,
            "clock_rate": corpus.truth.clock_rate,
        },
        "attestation_observation_types": _count_obs_types(corpus.attestations),
    }
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def _count_obs_types(attestations: list[Attestation]) -> dict[str, int]:
    out: dict[str, int] = {}
    for a in attestations:
        out[a.observation_type] = out.get(a.observation_type, 0) + 1
    return out
