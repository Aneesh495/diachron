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
    assignment_confidence: float
    quotation_excerpt: str | None
    quotation_truncated: bool
    conflicts_with: list[str]


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

        # Conflicting records: other attestations of the same language whose
        # dating intervals do not overlap this one's.
        conflicts = list(att.conflicting_with)
        for other in self._by_language.get(att.language_id, []):
            if other.id == att.id:
                continue
            if not self._intervals_overlap(
                (att.age_low, att.age_high), (other.age_low, other.age_high)
            ):
                if other.id not in conflicts:
                    conflicts.append(other.id)

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
        )

    def _limited_quotation(self, quotation: str | None) -> tuple[str | None, bool]:
        if not quotation:
            return None, False
        words = quotation.split()
        if len(words) <= _MAX_QUOTE_WORDS:
            return quotation, False
        return " ".join(words[:_MAX_QUOTE_WORDS]) + " [...]", True

    @staticmethod
    def _intervals_overlap(a: tuple[float, float], b: tuple[float, float]) -> bool:
        return not (a[1] < b[0] or b[1] < a[0])

    def conflicts_for_language(self, language_id: str) -> list[tuple[str, str]]:
        """Pairs of attestation ids for one language whose dates do not overlap."""
        records = self._by_language.get(language_id, [])
        out = []
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                if not self._intervals_overlap(
                    (records[i].age_low, records[i].age_high),
                    (records[j].age_low, records[j].age_high),
                ):
                    out.append((records[i].id, records[j].id))
        return out
