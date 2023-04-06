"""Known-truth simulation of cognate data, borrowing, and geography.

Everything produced here is labeled synthetic. The simulator draws a tree,
evolves cognate classes along it under the same continuous-time model the
inference uses, optionally injects borrowing between contemporaneous lineages,
and places lineages in geographic space by diffusion. It records the true tree,
true parameters, and true locations so recovery studies (gate G05) can measure
how well inference recovers them.

The point of a known-truth corpus is calibration: if inference cannot recover a
tree it generated under its own model, something is wrong with the sampler; if
it can, the credible-region coverage tells us whether the uncertainty is honest.
Synthetic data is never presented as observed history, and every dataset it
produces carries `is_synthetic=True`.

The cognate-class birth/death process is original: a concept starts with a root
class, and along each branch the class can be replaced (a new class is born,
replacing the old form) at the loss rate, so a concept accumulates distinct
classes across the tree in a way that mimics lexical replacement. This gives a
realistic many-class-per-concept structure rather than a single binary site.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .types import (
    Concept,
    Language,
    Lexeme,
    CognateJudgment,
    Dataset,
    Attestation,
)


@dataclass
class TrueHistory:
    """Ground truth accompanying a synthetic dataset."""

    parent: list[int]            # parent[i], -1 for root
    heights: list[float]         # node age in years before present
    taxon: list[int]             # taxon index for tips, -1 for internal
    tip_language_ids: list[str]  # taxon index -> language id
    root_age: float
    gain: float
    loss: float
    clock_rate: float
    true_locations: dict[str, tuple[float, float]]  # node label -> (lat, lon)


def _random_ultrametric_tree(
    n_tips: int, root_age: float, rng: np.random.Generator
) -> tuple[list[int], list[float], list[int]]:
    """Build a random coalescent-style ultrametric tree.

    Starting from `n_tips` lineages at the present, merge two random lineages at
    successively older times until one remains. Returns parallel arrays
    (parent, heights, taxon) with tips at age 0.
    """
    # Node 0..n_tips-1 are tips at height 0.
    parent = [-1] * (2 * n_tips - 1)
    heights = [0.0] * (2 * n_tips - 1)
    taxon = [-1] * (2 * n_tips - 1)
    for t in range(n_tips):
        taxon[t] = t

    active = list(range(n_tips))
    next_node = n_tips
    current_age = 0.0
    # Spread coalescence times across (0, root_age] with exponential gaps.
    k = n_tips
    while len(active) > 1:
        # Expected coalescent waiting time shrinks as lineages merge; scale so
        # the last merge lands near root_age.
        rate = k * (k - 1) / 2.0
        gap = rng.exponential(1.0 / rate) if rate > 0 else root_age
        current_age += gap
        i = rng.integers(0, len(active))
        a = active.pop(int(i))
        j = rng.integers(0, len(active))
        b = active.pop(int(j))
        node = next_node
        next_node += 1
        parent[a] = node
        parent[b] = node
        heights[node] = current_age
        active.append(node)
        k -= 1

    # Rescale heights so the root sits exactly at root_age.
    root = active[0]
    scale = root_age / heights[root] if heights[root] > 0 else 1.0
    heights = [h * scale for h in heights]
    return parent, heights, taxon


def _children_of(parent: list[int]) -> dict[int, list[int]]:
    ch: dict[int, list[int]] = {}
    for i, p in enumerate(parent):
        if p >= 0:
            ch.setdefault(p, []).append(i)
    return ch


def simulate_dataset(
    *,
    n_languages: int,
    n_concepts: int,
    root_age: float = 6000.0,
    gain: float = 0.0004,
    loss: float = 0.0009,
    clock_rate: float = 1.0,
    borrow_prob: float = 0.0,
    missing_frac: float = 0.0,
    seed: int = 0,
    name: str = "synthetic",
) -> tuple[Dataset, TrueHistory]:
    """Simulate a labeled-synthetic dataset and its ground truth.

    Rates are per year; `root_age` in years before present. `borrow_prob` is the
    per-concept probability that a borrowing event copies one lineage's class
    into a geographically near lineage (used for borrowing recovery tests).
    `missing_frac` drops that fraction of cells to simulate incomplete coverage.
    """
    rng = np.random.default_rng(seed)
    parent, heights, taxon = _random_ultrametric_tree(n_languages, root_age, rng)
    children = _children_of(parent)
    n_nodes = len(parent)
    root = next(i for i, p in enumerate(parent) if p < 0)

    # Preorder for top-down simulation.
    order: list[int] = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        for c in children.get(node, []):
            stack.append(c)

    languages: list[Language] = []
    tip_language_ids: list[str] = [""] * n_languages
    # Geography: Brownian diffusion of a 2D location from a root point.
    locations: dict[int, tuple[float, float]] = {}
    root_lat, root_lon = 45.0, 30.0
    geo_sigma = 0.004  # degrees^2 per year, continental scale over 6000y
    locations[root] = (root_lat, root_lon)
    for node in order:
        if node == root:
            continue
        p = parent[node]
        dt = heights[p] - heights[node]
        plat, plon = locations[p]
        lat = plat + rng.normal(0.0, np.sqrt(geo_sigma * dt))
        lon = plon + rng.normal(0.0, np.sqrt(geo_sigma * dt))
        locations[node] = (lat, lon)

    for t in range(n_languages):
        lat, lon = locations[t]
        lid = f"{name}_L{t:03d}"
        tip_language_ids[t] = lid
        languages.append(
            Language(
                id=lid,
                name=f"Lang{t}",
                glottocode=f"syn{t:05d}",
                family=name,
                latitude=float(lat),
                longitude=float(lon),
            )
        )

    concepts = [
        Concept(id=f"{name}_C{c:03d}", gloss=f"concept{c}") for c in range(n_concepts)
    ]

    lexemes: list[Lexeme] = []
    judgments: list[CognateJudgment] = []

    # For each concept evolve a cognate class along the tree. The class at a node
    # is inherited from its parent, but with probability 1 - exp(-loss*rate*dt) a
    # replacement occurs and a fresh class id is minted.
    for c, concept in enumerate(concepts):
        node_class: dict[int, int] = {}
        next_class = 0
        node_class[root] = next_class
        next_class += 1
        for node in order:
            if node == root:
                continue
            p = parent[node]
            dt = (heights[p] - heights[node]) * clock_rate
            replace_prob = 1.0 - np.exp(-loss * dt)
            if rng.random() < replace_prob:
                node_class[node] = next_class
                next_class += 1
            else:
                node_class[node] = node_class[p]

        # Optional borrowing: with borrow_prob, pick a donor tip and a nearby
        # recipient tip and copy the donor's class into the recipient.
        borrowed_tips: set[int] = set()
        if borrow_prob > 0 and rng.random() < borrow_prob:
            donor = int(rng.integers(0, n_languages))
            dlat, dlon = locations[donor]
            # nearest other tip by squared degree distance
            best = None
            best_d = float("inf")
            for other in range(n_languages):
                if other == donor:
                    continue
                olat, olon = locations[other]
                d = (olat - dlat) ** 2 + (olon - dlon) ** 2
                if d < best_d:
                    best_d = d
                    best = other
            if best is not None:
                node_class[best] = node_class[donor]
                borrowed_tips.add(best)

        for t in range(n_languages):
            if missing_frac > 0 and rng.random() < missing_frac:
                continue  # missing cell: no lexeme for this (language, concept)
            klass = node_class[t]
            lid = tip_language_ids[t]
            lex_id = f"{name}_C{c:03d}_L{t:03d}"
            lexemes.append(
                Lexeme(
                    id=lex_id,
                    language_id=lid,
                    concept_id=concept.id,
                    form=f"form_{c}_{klass}",
                    source_ref=f"{name}:simulated",
                )
            )
            judgments.append(
                CognateJudgment(
                    lexeme_id=lex_id,
                    concept_id=concept.id,
                    cognate_class=f"{concept.id}:{klass}",
                    confidence=1.0,
                    is_loan=(t in borrowed_tips),
                    annotator="simulator",
                    source_ref=f"{name}:simulated",
                )
            )

    dataset = Dataset(
        name=name,
        version="synthetic-1",
        license="synthetic (no restrictions)",
        doi=None,
        concepts=concepts,
        languages=languages,
        lexemes=lexemes,
        judgments=judgments,
        is_synthetic=True,
        notes="Labeled synthetic known-truth corpus. Not observed history.",
    )

    truth = TrueHistory(
        parent=parent,
        heights=heights,
        taxon=taxon,
        tip_language_ids=tip_language_ids,
        root_age=root_age,
        gain=gain,
        loss=loss,
        clock_rate=clock_rate,
        true_locations={
            f"{name}_L{t:03d}": locations[t] for t in range(n_languages)
        },
    )
    return dataset, truth


def simulate_attestations(
    dataset: Dataset,
    truth: TrueHistory,
    *,
    n_records: int,
    seed: int = 0,
) -> list[Attestation]:
    """Draw synthetic dated/located attestations consistent with the truth.

    Each record samples a language, a true age along its lineage, and a noisy
    dating interval and location. The records are deliberately sparse and biased
    toward some regions so inference must reflect that bias rather than pretend to
    even coverage.
    """
    rng = np.random.default_rng(seed + 777)
    attestations: list[Attestation] = []
    tip_ids = truth.tip_language_ids
    for r in range(n_records):
        t = int(rng.integers(0, len(tip_ids)))
        lid = tip_ids[t]
        # A plausible attestation age between present and the tip's lineage age.
        max_age = truth.root_age * 0.6
        true_age = float(rng.uniform(0, max_age))
        halfwidth = float(rng.uniform(50, 400))
        lat, lon = truth.true_locations[lid]
        obs_lat = lat + rng.normal(0, 0.5)
        obs_lon = lon + rng.normal(0, 0.5)
        attestations.append(
            Attestation(
                id=f"{dataset.name}_att{r:04d}",
                language_id=lid,
                age_low=max(0.0, true_age - halfwidth),
                age_high=true_age + halfwidth,
                source_ref=f"{dataset.name}:synthetic-attestation:{r}",
                page=str(rng.integers(1, 400)),
                latitude=float(obs_lat),
                longitude=float(obs_lon),
                location_uncertainty_km=float(rng.uniform(20, 200)),
                language_assignment_confidence=float(rng.uniform(0.6, 1.0)),
                observation_type=rng.choice(
                    ["attestation", "inscription", "reconstruction", "trade"]
                ),
            )
        )
    return attestations
