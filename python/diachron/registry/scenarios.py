"""Alternative-history registry.

Every historical scenario Diachron entertains is recorded here against the exact
model, assumptions, evidence, and the posterior mass or sampling frequency that
supports it. A scenario is a falsifiable claim ("the family's root is older than
7000 years", "clade {A,B,C} is monophyletic", "the homeland lies in region R"),
and the registry attaches to it the model configuration under which it was
evaluated and the posterior support measured from saved draws.

The point is that no scenario floats free of its evidence: the registry refuses
to store a scenario without a model identity and a support value computed from a
posterior bundle, so a claim can always be traced to what produced it and to the
assumptions it depends on.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
import json
import math
import os
from typing import Callable

from ..infer.artifacts import PosteriorBundle
from ..diagnostics import topology as tp
from ..newick_support import clade_support_for  # local helper, see below


@dataclass
class Scenario:
    """A falsifiable historical claim and the support measured for it."""

    name: str
    description: str
    model_name: str
    assumptions: list[str]
    evidence_refs: list[str]
    support: float               # posterior probability or sampling frequency
    support_kind: str            # "clade_frequency" | "scalar_interval" | "region_mass"
    detail: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


class ScenarioRegistry:
    """Collects scenarios and writes them with their evidence."""

    def __init__(self) -> None:
        self.scenarios: list[Scenario] = []

    def add(self, scenario: Scenario) -> None:
        if not scenario.model_name:
            raise ValueError("a scenario must name the model it was evaluated under")
        if not math.isfinite(scenario.support) or not 0.0 <= scenario.support <= 1.0:
            raise ValueError("scenario support must be a finite sampling frequency")
        if scenario.support_kind not in {
            "clade_frequency",
            "scalar_interval",
            "region_mass",
        }:
            raise ValueError(f"unknown support kind {scenario.support_kind}")
        self.scenarios.append(scenario)

    # --- scenario builders tied to a posterior bundle --------------------
    def add_clade_scenario(
        self,
        bundle: PosteriorBundle,
        taxa: set[int],
        name: str,
        description: str,
        assumptions: list[str],
        evidence_refs: list[str],
    ) -> Scenario:
        """A monophyly claim: posterior frequency of a clade over all draws."""
        newicks = [
            nw for c in range(bundle.n_chains) for nw in bundle.chain_newicks(c)
        ]
        posteriors = tp.clade_posteriors(newicks)
        if not newicks or len(taxa) < 2 or any(not 0 <= t < len(bundle.tip_labels) for t in taxa):
            raise ValueError("clade scenario requires draws and at least two known taxa")
        labels = {bundle.tip_labels[t] for t in taxa}
        support = clade_support_for(posteriors, labels)
        sc = Scenario(
            name=name,
            description=description,
            model_name=bundle.model_name,
            assumptions=assumptions,
            evidence_refs=evidence_refs,
            support=support,
            support_kind="clade_frequency",
            detail={"taxa": sorted(labels), **_identity(bundle)},
        )
        self.add(sc)
        return sc

    def add_age_scenario(
        self,
        bundle: PosteriorBundle,
        name: str,
        description: str,
        threshold_years: float,
        direction: str,
        assumptions: list[str],
        evidence_refs: list[str],
    ) -> Scenario:
        """A chronology claim: posterior probability the root age exceeds/precedes a threshold."""
        ages = [d["root_age"] for chain in bundle.chains for d in chain]
        if not ages or not math.isfinite(threshold_years) or threshold_years < 0:
            raise ValueError("age scenario requires draws and a finite nonnegative age")
        if direction == "older_than":
            support = sum(1 for a in ages if a > threshold_years) / len(ages)
        elif direction == "younger_than":
            support = sum(1 for a in ages if a < threshold_years) / len(ages)
        else:
            raise ValueError("direction must be older_than or younger_than")
        sc = Scenario(
            name=name,
            description=description,
            model_name=bundle.model_name,
            assumptions=assumptions,
            evidence_refs=evidence_refs,
            support=support,
            support_kind="scalar_interval",
            detail={"threshold_years": threshold_years, "direction": direction, **_identity(bundle)},
        )
        self.add(sc)
        return sc

    def add_custom_scenario(
        self,
        bundle: PosteriorBundle,
        name: str,
        description: str,
        support_fn: Callable[[PosteriorBundle], float],
        support_kind: str,
        assumptions: list[str],
        evidence_refs: list[str],
    ) -> Scenario:
        sc = Scenario(
            name=name,
            description=description,
            model_name=bundle.model_name,
            assumptions=assumptions,
            evidence_refs=evidence_refs,
            support=float(support_fn(bundle)),
            support_kind=support_kind,
            detail=_identity(bundle),
        )
        self.add(sc)
        return sc

    def save(self, path: str) -> None:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump([s.to_dict() for s in self.scenarios], f, indent=2)

    def summary(self) -> list[dict]:
        return [
            {
                "name": s.name,
                "model": s.model_name,
                "support": round(s.support, 3),
                "kind": s.support_kind,
            }
            for s in self.scenarios
        ]


def _identity(bundle: PosteriorBundle) -> dict:
    return {
        "dataset": bundle.dataset_name, "dataset_version": bundle.dataset_version,
        "is_synthetic": bundle.is_synthetic, "source_commit": bundle.source_commit,
        "draw_digest": bundle.draw_digest(), "n_draws": sum(len(c) for c in bundle.chains),
        "config": bundle.context.get("model", {}).get("config"),
        "scope": "Sampling frequency in the named bundle.",
    }
