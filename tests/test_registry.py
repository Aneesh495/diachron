"""Tests for the alternative-history registry and attestation inspector."""

from __future__ import annotations

import pytest

from diachron.io.types import Attestation
from diachron.registry import ScenarioRegistry, Scenario, AttestationInspector


def test_registry_refuses_scenario_without_model():
    reg = ScenarioRegistry()
    with pytest.raises(ValueError):
        reg.add(
            Scenario(
                name="x",
                description="d",
                model_name="",
                assumptions=[],
                evidence_refs=[],
                support=0.5,
                support_kind="clade_frequency",
            )
        )


def test_registry_rejects_unknown_support_kind():
    reg = ScenarioRegistry()
    with pytest.raises(ValueError):
        reg.add(
            Scenario(
                name="x",
                description="d",
                model_name="dated",
                assumptions=[],
                evidence_refs=[],
                support=0.5,
                support_kind="made_up",
            )
        )


def test_registry_saves_and_summarizes(tmp_path):
    reg = ScenarioRegistry()
    reg.add(
        Scenario(
            name="root_old",
            description="root older than 6000",
            model_name="dated",
            assumptions=["dated model"],
            evidence_refs=["manifest"],
            support=0.9,
            support_kind="scalar_interval",
        )
    )
    path = str(tmp_path / "scenarios.json")
    reg.save(path)
    summ = reg.summary()
    assert summ[0]["model"] == "dated"
    assert summ[0]["support"] == 0.9


def _atts():
    return [
        Attestation(
            id="a1", language_id="L1", age_low=100, age_high=300,
            source_ref="ref1", quotation=" ".join(["word"] * 40), conflicting_with=("a2",),
        ),
        # conflicting with a1 (non-overlapping interval) for the same language
        Attestation(
            id="a2", language_id="L1", age_low=500, age_high=700, source_ref="ref2",
        ),
        Attestation(
            id="a4", language_id="L1", age_low=900, age_high=950, source_ref="independent event",
        ),
        Attestation(
            id="a3", language_id="L2", age_low=100, age_high=200, source_ref="ref3",
        ),
    ]


def test_inspector_truncates_quotation():
    insp = AttestationInspector(_atts())
    rep = insp.inspect("a1")
    assert rep.quotation_truncated
    assert rep.quotation_excerpt.endswith("[...]")
    assert len(rep.quotation_excerpt.split()) <= 26  # 25 words + marker


def test_inspector_surfaces_conflicts():
    insp = AttestationInspector(_atts())
    rep = insp.inspect("a1")
    assert "a2" in rep.conflicts_with
    pairs = insp.conflicts_for_language("L1")
    assert ("a1", "a2") in pairs


def test_inspector_no_false_conflict():
    insp = AttestationInspector(_atts())
    # L2 has a single attestation; no conflicts.
    assert insp.conflicts_for_language("L2") == []
    assert insp.inspect("a4").conflicts_with == []


def test_scenario_uses_actual_bundle_taxon_ids():
    from diachron.infer.artifacts import PosteriorBundle

    bundle = PosteriorBundle(
        "dated", "observed", "1", False, 1, 1, ["var-A", "var-B", "var-C"],
        [[{"newick": "((var-A:1,var-B:1):1,var-C:2);"}]], source_commit="example",
    )
    scenario = ScenarioRegistry().add_clade_scenario(
        bundle, {0, 1}, "sister-pair", "A and B form a clade", ["sampled tree"], ["bundle"],
    )
    assert scenario.support == 1.0
    assert scenario.detail["taxa"] == ["var-A", "var-B"]
    assert scenario.detail["source_commit"] == "example"
