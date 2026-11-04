"""Multi-chain sampling campaign and the data-to-model wiring.

A campaign runs several independently initialized chains for one model on one
dataset, then packages their draws into a posterior bundle. Independent
initialization means each chain starts from its own random tree and its own RNG
seed, so between-chain agreement is real evidence of convergence rather than an
artifact of a shared start.

The module also builds the `PosteriorModel` from a dataset: it transforms cognate
judgments to characters, aligns the character columns to the tree's taxon order,
derives node-age calibrations from dated attestations, and attaches tip
coordinates for the geographic model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from ..io.types import Dataset, Attestation
from ..io.characters import build_character_matrix
from ..models.tree import random_starting_tree
from ..models.tip_dates import TipDate
from ..models.posterior import PosteriorModel, ModelConfig, Calibration
from .chain import Chain, ChainConfig, ChainState
from .artifacts import bundle_from_chains, PosteriorBundle


@dataclass
class StudyData:
    """Everything a campaign needs derived from a dataset."""

    tip_labels: list[str]          # taxon index -> language id
    columns: list[list[int]]       # variable character columns in taxon order
    calibrations: list[Calibration]
    tip_coords: dict[int, tuple[float, float]]
    dataset: Dataset
    borrow_candidates: list[tuple[int, int]] = None  # type: ignore[assignment]
    minimum_present: int = 1
    attestations: list[Attestation] = field(default_factory=list)
    tip_dates: dict[int, TipDate] = field(default_factory=dict)
    character_ids: list[tuple[str, str]] = field(default_factory=list)


def prepare_study(
    dataset: Dataset,
    attestations: list[Attestation] | None = None,
    *,
    drop_singletons: bool = True,
    include_loans_as_present: bool = True,
    include_unobserved_languages: bool = False,
) -> StudyData:
    """Transform a dataset into model-ready arrays aligned to a taxon order."""
    cm = build_character_matrix(
        dataset,
        drop_singletons=drop_singletons,
        include_loans_as_present=include_loans_as_present,
        include_unobserved_languages=include_unobserved_languages,
    )
    tip_labels = list(cm.language_ids)
    lang_to_idx = {lid: i for i, lid in enumerate(tip_labels)}

    characters = cm.variable_characters()
    columns = [char.states for char in characters]

    # Tip coordinates for the geographic model, keyed by taxon index.
    coord_by_lang = {
        lang.id: (lang.latitude, lang.longitude)
        for lang in dataset.languages if not lang.is_historical
    }
    tip_coords: dict[int, tuple[float, float]] = {}
    for t, lid in enumerate(tip_labels):
        lat, lon = coord_by_lang.get(lid, (None, None))
        if lat is not None and lon is not None:
            tip_coords[t] = (float(lat), float(lon))

    # Calibrations from attestations. A single contemporary tip gives no
    # time-scale information, so we do not place one calibration per attestation.
    # Instead the attestations collectively bound the age of the whole family:
    # the oldest reliable attestation is a lower bound on the root age, and a
    # defended upper bound (here taken a generation beyond the oldest, reflecting
    # that the proto-language predates its oldest attested descendant) caps it.
    # This yields one identifiable root-age calibration (the MRCA of all tips),
    # which is the minimal structure that makes the strict clock's time scale
    # identifiable without over-constraining every internal node.
    calibrations: list[Calibration] = []
    usable = [
        a for a in (attestations or [])
        if a.language_id in lang_to_idx
        and a.language_assignment_confidence is not None
        and a.language_assignment_confidence >= 0.5
    ]
    if usable and dataset.is_synthetic:
        all_taxa = tuple(range(len(tip_labels)))
        oldest = max(a.age_high for a in usable)
        # Lower bound: the oldest attested age. Upper bound: a margin above it,
        # because the common ancestor must predate the oldest attested descendant.
        root_low = oldest
        root_high = oldest * 1.8 + 500.0
        calibrations.append(
            Calibration(
                clade_taxa=all_taxa,
                age_low=root_low,
                age_high=root_high,
                kind="uniform",
                source_ref="aggregate of dated attestations",
            )
        )

    # Bounded borrowing candidates: each tip's nearest geographic neighbor, as a
    # fixed, finite candidate set. This is the "bounded network" the borrowing
    # model operates on; it is derived from modern reference coordinates and is
    # explicitly a candidate set, not a claim of historical contact.
    borrow_candidates = _nearest_neighbor_candidates(tip_coords, max_edges=20)

    tip_dates = {
        lang_to_idx[l.id]: TipDate(l.date_distribution, l.date_mean, l.date_sd, l.date_offset)
        for l in dataset.languages if l.id in lang_to_idx and l.date_distribution is not None
    }
    return StudyData(
        tip_labels=tip_labels,
        columns=columns,
        calibrations=calibrations,
        tip_coords=tip_coords,
        dataset=dataset,
        borrow_candidates=borrow_candidates,
        minimum_present=2 if drop_singletons else 1,
        attestations=list(attestations or ()),
        tip_dates=tip_dates,
        character_ids=[(char.concept_id, char.cognate_class) for char in characters],
    )


def _nearest_neighbor_candidates(
    tip_coords: dict[int, tuple[float, float]], max_edges: int
) -> list[tuple[int, int]]:
    """Each tip paired with its nearest geographic neighbor, capped to max_edges.

    Uses squared-degree distance (cheap and monotone with great-circle at these
    scales) and deduplicates symmetric pairs so (a,b) and (b,a) count once.
    """
    taxa = sorted(tip_coords)
    seen: set[tuple[int, int]] = set()
    edges: list[tuple[int, int, float]] = []
    for a in taxa:
        lat_a, lon_a = tip_coords[a]
        best = None
        best_d = float("inf")
        for b in taxa:
            if b == a:
                continue
            lat_b, lon_b = tip_coords[b]
            d = (lat_a - lat_b) ** 2 + (lon_a - lon_b) ** 2
            if d < best_d:
                best_d = d
                best = b
        if best is not None:
            key = (min(a, best), max(a, best))
            if key not in seen:
                seen.add(key)
                edges.append((key[0], key[1], best_d))
    edges.sort(key=lambda e: e[2])
    return [(a, b) for a, b, _ in edges[:max_edges]]


def build_model(study: StudyData, config: ModelConfig) -> PosteriorModel:
    """Assemble a PosteriorModel for one configuration from prepared study data."""
    # Calibrations and coordinates are only attached when the config uses them,
    # but the model stores them regardless and the config flags gate their use.
    candidates = study.borrow_candidates or []
    return PosteriorModel(
        config=config,
        columns=study.columns,
        n_states=2,
        minimum_present=study.minimum_present,
        tip_dates=study.tip_dates if config.use_calibration else {},
        calibrations=study.calibrations,
        tip_coords=study.tip_coords,
        borrow_candidates=candidates if config.use_borrowing else [],
    )


def run_campaign(
    study: StudyData,
    config: ModelConfig,
    *,
    n_chains: int = 4,
    warmup: int = 1500,
    iterations: int = 6000,
    thin: int = 10,
    root_age_init: float = 6000.0,
    base_seed: int = 1000,
    checkpoint_dir: str | None = None,
    progress=None,
    resume: bool = False,
) -> list[Chain]:
    """Run `n_chains` independently initialized chains for one model."""
    model = build_model(study, config)
    chains: list[Chain] = []
    for c in range(n_chains):
        rng = np.random.default_rng(base_seed + c)
        dates = {taxon: date.center() for taxon, date in model.tip_dates.items()}
        start = random_starting_tree(study.tip_labels, root_age_init, rng, dates)
        cfg = ChainConfig(
            seed=base_seed + 7919 * (c + 1),
            warmup=warmup,
            iterations=iterations,
            thin=thin,
        )
        borrow_init = (
            [0.02] * len(study.borrow_candidates)
            if config.use_borrowing and study.borrow_candidates
            else []
        )
        state = ChainState(
            gain=0.0005, loss=0.001, clock_rate=1.0, borrow_rates=borrow_init
        )
        ckpt = (
            Path(checkpoint_dir) / f"{config.name}_chain{c}.json" if checkpoint_dir else None
        )
        if ckpt is not None and ckpt.exists():
            if not resume:
                raise FileExistsError(f"checkpoint already exists: {ckpt}; use --resume")
            chain = Chain.resume(model, str(ckpt))
            if chain.config != cfg or chain.tree.tip_labels != start.tip_labels:
                raise ValueError(f"checkpoint chain budget, seed, or tip labels differ: {ckpt}")
        else:
            chain = Chain(model, start, state, cfg)
        chain.run(checkpoint_path=ckpt, checkpoint_every=1000)
        chains.append(chain)
        if progress:
            progress(config.name, c, chain)
    return chains


def campaign_bundle(
    chains: list[Chain],
    study: StudyData,
    config: ModelConfig,
    source_commit: str = "unknown",
) -> PosteriorBundle:
    bundle = bundle_from_chains(
        chains,
        model_name=config.name,
        dataset_name=study.dataset.name,
        dataset_version=study.dataset.version,
        is_synthetic=study.dataset.is_synthetic,
        tip_labels=study.tip_labels,
        source_commit=source_commit,
    )
    bundle.context = {
        "schema_version": 2,
        "model": asdict(chains[0].model),
        "attestations": [asdict(a) for a in study.attestations],
        "dataset_url": study.dataset.source_url,
    }
    return bundle
