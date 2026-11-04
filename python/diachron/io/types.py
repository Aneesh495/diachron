"""Typed records for lexical data, cognate judgments, and attestations.

These dataclasses are the in-memory contract between ingestion and the models.
They are deliberately explicit about the distinctions the domain brief insists
on: a missing observation is not an absence, a loan is annotated rather than
silently treated as a cognate, and multiple conflicting judgments for one cell
are preserved rather than collapsed.

Identifiers are strings carried verbatim from the source (Glottolog codes,
CLDF ids, concept ids). They survive renamed or merged catalog entries because a
frozen dataset records the identifiers it was built with; downstream code never
rewrites them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


class DataError(Exception):
    """Raised when ingested data violates a stated contract."""


@dataclass(frozen=True)
class Concept:
    """A meaning slot, e.g. the Swadesh concept 'water'.

    `concepticon_id` links to a cross-dataset concept catalog when available so
    two datasets can be aligned on meaning rather than on a gloss string.
    """

    id: str
    gloss: str
    concepticon_id: Optional[str] = None
    semantic_field: Optional[str] = None


@dataclass(frozen=True)
class Language:
    """A language or dialect variety.

    `glottocode` is the versioned Glottolog identifier. `latitude`/`longitude`
    are the catalog reference point, explicitly a modern reference location and
    not a historical homeland; geographic inference treats them as uncertain
    observations, never as boundaries.
    """

    id: str
    name: str
    glottocode: Optional[str] = None
    family: Optional[str] = None
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    is_extinct: bool = False
    # Age of the variety's attestation in years before present, when the tip is
    # dated (e.g. an extinct documented language). None for contemporary tips.
    tip_age: Optional[float] = None
    is_historical: bool = False
    source_classification: Optional[str] = None
    date_distribution: Optional[str] = None
    date_mean: Optional[float] = None
    date_sd: Optional[float] = None
    date_offset: Optional[float] = None


@dataclass(frozen=True)
class Lexeme:
    """A single observed form for a (language, concept) cell.

    `form` is the transliterated form used for display and alignment; `original`
    preserves the source script. `transliteration_source` records how `form` was
    derived so orthographic provenance is never lost. Orthographic similarity is
    never equated with cognacy; that is what the cognate judgments are for.
    """

    id: str
    language_id: str
    concept_id: str
    form: str
    original: Optional[str] = None
    transliteration_source: Optional[str] = None
    source_ref: Optional[str] = None
    source_value: Optional[str] = None
    loan_annotation: Optional[bool] = None


@dataclass(frozen=True)
class CognateJudgment:
    """An assignment of a lexeme to a cognate class for a concept.

    `cognate_class` groups lexemes judged to descend from a common ancestral
    form for that concept. `confidence` in [0, 1] expresses uncertainty in the
    judgment; `is_loan` marks a form believed to be borrowed rather than
    inherited; `is_partial` marks a partial cognate (shared morpheme). Multiple
    judgments may disagree for one lexeme, and all are retained with their
    `annotator` so conflicts are visible.
    """

    lexeme_id: str
    concept_id: str
    cognate_class: str
    confidence: float = 1.0
    is_loan: bool = False
    is_partial: bool = False
    annotator: Optional[str] = None
    source_ref: Optional[str] = None
    id: Optional[str] = None
    doubtful: Optional[bool] = None
    segment_slice: Optional[str] = None
    detection_method: Optional[str] = None
    loan_event: Optional[str] = None
    loan_source: Optional[str] = None


@dataclass(frozen=True)
class CognateClassRecord:
    """A source class and optional alternative links, before binary projection."""

    id: str
    root_form: str | None
    root_language: str | None
    source_ref: str | None
    proposed_cognate_id: str | None
    proposed_scale: int | None
    superset_id: str | None
    parallel_derivation: bool | None


@dataclass(frozen=True)
class Attestation:
    """A dated, located record of a language at a point in the past.

    An attestation is an observation with its own uncertainty: a dated inscription,
    a documented text, a reconstructed form with a defended date. It is never a
    whole-population distribution. `age_low`/`age_high` bracket the date in years
    before present (a conflicting pair of records produces two attestations, not a
    merged midpoint). `region` and `location_uncertainty_km` describe where, with
    the understanding that an elite inscription or a trade document need not mark
    the centroid of a speaker population. `observation_type` records which of
    these it is so the geographic likelihood can weight it appropriately.
    """

    id: str
    language_id: str
    age_low: float
    age_high: float
    source_ref: str
    page: Optional[str] = None
    region: Optional[str] = None
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    location_uncertainty_km: Optional[float] = None
    language_assignment_confidence: Optional[float] = 1.0
    observation_type: str = "attestation"  # attestation | inscription | reconstruction | trade
    conflicting_with: tuple[str, ...] = field(default_factory=tuple)
    quotation: Optional[str] = None  # exposed only within permitted limits
    bibliography: tuple[str, ...] = ()
    location_kind: Optional[str] = None
    location_uncertainty_note: Optional[str] = None
    language_assignment_note: Optional[str] = None
    source_date_start: Optional[str] = None
    source_date_end: Optional[str] = None


@dataclass
class Dataset:
    """A frozen bundle of lexical data for one study.

    The dataset records the source identity (name, version, DOI, license) so a
    frozen study can be audited against the exact release it was built from.
    `is_synthetic` is True for labeled known-truth corpora and is propagated into
    every artifact and view so synthetic data is never presented as observed.
    """

    name: str
    version: str
    license: str
    doi: Optional[str]
    concepts: list[Concept]
    languages: list[Language]
    lexemes: list[Lexeme]
    judgments: list[CognateJudgment]
    is_synthetic: bool = False
    source_url: Optional[str] = None
    notes: str = ""
    cognate_classes: list[CognateClassRecord] = field(default_factory=list)

    def validate(self) -> None:
        """Check referential integrity and the invariants the models rely on."""
        concept_ids = {c.id for c in self.concepts}
        language_ids = {lang.id for lang in self.languages}
        lexeme_ids = {lx.id for lx in self.lexemes}

        if len(concept_ids) != len(self.concepts):
            raise DataError("duplicate concept ids")
        if len(language_ids) != len(self.languages):
            raise DataError("duplicate language ids")
        if len(lexeme_ids) != len(self.lexemes):
            raise DataError("duplicate lexeme ids")
        class_ids = {row.id for row in self.cognate_classes}
        if len(class_ids) != len(self.cognate_classes):
            raise DataError("duplicate source cognate class ids")

        for lx in self.lexemes:
            if lx.language_id not in language_ids:
                raise DataError(f"lexeme {lx.id} references unknown language {lx.language_id}")
            if lx.concept_id not in concept_ids:
                raise DataError(f"lexeme {lx.id} references unknown concept {lx.concept_id}")

        for j in self.judgments:
            if j.lexeme_id not in lexeme_ids:
                raise DataError(f"judgment references unknown lexeme {j.lexeme_id}")
            if j.concept_id not in concept_ids:
                raise DataError(f"judgment references unknown concept {j.concept_id}")
            if not (0.0 <= j.confidence <= 1.0):
                raise DataError(f"judgment confidence out of range: {j.confidence}")

    @property
    def n_eligible_languages(self) -> int:
        """Languages with at least one lexeme (an observed dimension)."""
        with_data = {lx.language_id for lx in self.lexemes}
        return sum(1 for lang in self.languages if lang.id in with_data)

    @property
    def n_eligible_concepts(self) -> int:
        with_data = {lx.concept_id for lx in self.lexemes}
        return sum(1 for c in self.concepts if c.id in with_data)

    def observed_cell_count(self) -> int:
        """Number of (language, concept) cells with at least one lexeme."""
        cells = {(lx.language_id, lx.concept_id) for lx in self.lexemes}
        return len(cells)
