"""Attestation inspector.

Exposes a dated attestation's bibliographic reference, dating interval and its
uncertainty, geographic observation type, language-assignment confidence, and any
conflicting records, while limiting any source quotation to a permitted length.
The quotation cap is deliberate: Diachron attributes sources and shows enough to
be checkable without reproducing copyrighted text, so the inspector truncates a
quotation to a short excerpt and records that it was truncated.

Conflicting records are surfaced rather than merged: two attestations that date
the same language differently are both shown, with their disagreement made
explicit, so a dating conflict is visible to the researcher instead of being
averaged away.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..io.types import Attestation

# Maximum words of a source quotation exposed at once, to respect licensing.
_MAX_QUOTE_WORDS = 25


@dataclass
class AttestationReport:
    id: str
    language_id: str
    age_interval_ybp: tuple[float, float]
    age_midpoint: float
    age_halfwidth: float
    observation_type: str
    source_ref: str
    page: str | None
    location: tuple[float, float] | None
    location_uncertainty_km: float | None
    assignment_confidence: float | None
    quotation_excerpt: str | None
    quotation_truncated: bool
    conflicts_with: list[str]
    bibliography: tuple[str, ...]
    region: str | None
    location_kind: str | None
    location_note: str | None
    assignment_note: str | None
    source_dates: tuple[str | None, str | None]


class AttestationInspector:
    def __init__(self, attestations: list[Attestation]) -> None:
        self._by_id = {a.id: a for a in attestations}
        self._by_language: dict[str, list[Attestation]] = {}
        for a in attestations:
            self._by_language.setdefault(a.language_id, []).append(a)

    def inspect(self, attestation_id: str) -> AttestationReport:
        att = self._by_id[attestation_id]
        midpoint = 0.5 * (att.age_low + att.age_high)
        halfwidth = 0.5 * (att.age_high - att.age_low)

        excerpt, truncated = self._limited_quotation(att.quotation)

        # Conflict links refer to alternative accounts of the same evidence.
        # Distinct inscriptions in different years can both be correct.
        conflicts = sorted({
            other.id for other in self._by_language.get(att.language_id, ())
            if other.id != att.id and (
                other.id in att.conflicting_with or att.id in other.conflicting_with
            )
        })

        location = None
        if att.latitude is not None and att.longitude is not None:
            location = (att.latitude, att.longitude)

        return AttestationReport(
            id=att.id,
            language_id=att.language_id,
            age_interval_ybp=(att.age_low, att.age_high),
            age_midpoint=midpoint,
            age_halfwidth=halfwidth,
            observation_type=att.observation_type,
            source_ref=att.source_ref,
            page=att.page,
            location=location,
            location_uncertainty_km=att.location_uncertainty_km,
            assignment_confidence=att.language_assignment_confidence,
            quotation_excerpt=excerpt,
            quotation_truncated=truncated,
            conflicts_with=conflicts,
            bibliography=att.bibliography, region=att.region,
            location_kind=att.location_kind, location_note=att.location_uncertainty_note,
            assignment_note=att.language_assignment_note,
            source_dates=(att.source_date_start, att.source_date_end),
        )

    def _limited_quotation(self, quotation: str | None) -> tuple[str | None, bool]:
        if not quotation:
            return None, False
        words = quotation.split()
        if len(words) <= _MAX_QUOTE_WORDS:
            return quotation, False
        return " ".join(words[:_MAX_QUOTE_WORDS]) + " [...]", True

    def conflicts_for_language(self, language_id: str) -> list[tuple[str, str]]:
        """Explicit source alternatives for attestations of the same language."""
        return sorted({
            tuple(sorted((record.id, other)))
            for record in self._by_language.get(language_id, ())
            for other in self.inspect(record.id).conflicts_with
        })
