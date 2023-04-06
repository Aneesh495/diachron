"""Transform cognate judgments into the character matrix the likelihood uses.

The unit of inheritance is a cognate class within a concept. For each concept we
enumerate its observed cognate classes; each class becomes one binary character
whose state for a language is:

* present (1): the language has a lexeme judged to belong to this class;
* absent (0): the language has lexical data for this concept but no form in this
  class;
* missing (-1): the language has no lexical data for this concept at all, so we
  cannot say whether the class is present or absent.

The present/absent/missing distinction is the crux of gate G04: a language with
data for a concept but not in this class contributes an informative 0, while a
language with no data for the concept contributes a non-informative missing that
the pruning likelihood marginalizes over.

Eligibility and ascertainment (D06): singleton classes (present in a single
language) carry no phylogenetic signal about internal structure and are excluded
from the likelihood but recorded in the evidence ledger. The ascertainment
correction later conditions on the fact that we only keep variable characters.

Loans and partial cognates (D02, D03): a loan judgment is not folded into the
inherited character. A borrowed form is recorded in the character's provenance
and, in the borrowing-aware model, enters a separate contact channel. Partial
cognates are kept as their own classes rather than merged into a full cognate,
so shared morphology is not overclaimed as full cognacy.

Every character carries an `EvidenceLink` back to the lexemes and judgments that
produced each cell, which is what makes the transform auditable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .types import Dataset, DataError

# State codes. These match the convention used by the Rust binding, where a
# negative state means missing.
PRESENT = 1
ABSENT = 0
MISSING = -1


@dataclass
class Character:
    """One binary presence/absence character for a cognate class."""

    concept_id: str
    cognate_class: str
    # states[language_index] in {PRESENT, ABSENT, MISSING}, aligned to the
    # CharacterMatrix language order.
    states: list[int]
    is_partial: bool = False
    # Language ids whose present state came from a form annotated as a loan.
    loan_languages: tuple[str, ...] = field(default_factory=tuple)

    def n_present(self) -> int:
        return sum(1 for s in self.states if s == PRESENT)

    def n_observed(self) -> int:
        return sum(1 for s in self.states if s != MISSING)

    def is_variable(self) -> bool:
        """A character is variable if at least two observed languages differ.

        Constant characters (all observed present, or all observed absent) carry
        no signal for tree structure and are what the ascertainment correction
        conditions out.
        """
        observed = [s for s in self.states if s != MISSING]
        if len(observed) < 2:
            return False
        return len(set(observed)) > 1


@dataclass
class EvidenceLink:
    """Audit trail from a character cell back to its source judgments."""

    concept_id: str
    cognate_class: str
    language_id: str
    state: int
    lexeme_ids: tuple[str, ...]
    judgment_sources: tuple[str, ...]


@dataclass
class CharacterMatrix:
    """The full set of characters with the language order they are aligned to."""

    language_ids: list[str]
    characters: list[Character]
    evidence: list[EvidenceLink]
    excluded_singletons: int = 0
    excluded_constant: int = 0

    def language_index(self) -> dict[str, int]:
        return {lid: i for i, lid in enumerate(self.language_ids)}

    def variable_characters(self) -> list[Character]:
        return [c for c in self.characters if c.is_variable()]

    def column_for_core(self, character: Character) -> list[int]:
        """Return the per-language state list in the Rust binding's convention.

        The binding already treats negative states as missing, and PRESENT/ABSENT
        are 1/0, so the states list is passed through unchanged.
        """
        return list(character.states)


def build_character_matrix(
    dataset: Dataset,
    *,
    drop_singletons: bool = True,
    include_loans_as_present: bool = True,
) -> CharacterMatrix:
    """Transform a dataset into a character matrix with an evidence trail.

    `drop_singletons` excludes classes present in a single language from the
    returned characters (they are counted in `excluded_singletons`).
    `include_loans_as_present` controls whether a loan-annotated form sets the
    inherited character to present; the borrowing sensitivity analysis flips this
    to measure the effect of removing loans.
    """
    dataset.validate()

    languages = sorted({lx.language_id for lx in dataset.lexemes})
    lang_index = {lid: i for i, lid in enumerate(languages)}
    n_lang = len(languages)

    # Which (language, concept) cells have any lexical data? A language with no
    # lexeme for a concept is missing for every class of that concept.
    concept_langs: dict[str, set[str]] = {}
    for lx in dataset.lexemes:
        concept_langs.setdefault(lx.concept_id, set()).add(lx.language_id)

    # Index judgments by (concept, cognate_class) -> {language -> [lexeme ids]}.
    # Only judgments whose lexeme and concept agree are used.
    lexeme_by_id = {lx.id: lx for lx in dataset.lexemes}
    class_members: dict[tuple[str, str], dict[str, list[str]]] = {}
    class_is_partial: dict[tuple[str, str], bool] = {}
    class_loans: dict[tuple[str, str], set[str]] = {}
    judgment_src: dict[tuple[str, str, str], list[str]] = {}

    for j in dataset.judgments:
        lx = lexeme_by_id.get(j.lexeme_id)
        if lx is None or lx.concept_id != j.concept_id:
            raise DataError(f"judgment {j.lexeme_id} inconsistent with its lexeme")
        if j.is_loan and not include_loans_as_present:
            # Loan removal mode: the form does not count toward the inherited
            # class, but the cell is still observed (the language has data).
            continue
        key = (j.concept_id, j.cognate_class)
        members = class_members.setdefault(key, {})
        members.setdefault(lx.language_id, []).append(lx.id)
        class_is_partial[key] = class_is_partial.get(key, False) or j.is_partial
        if j.is_loan:
            class_loans.setdefault(key, set()).add(lx.language_id)
        src_key = (j.concept_id, j.cognate_class, lx.language_id)
        if j.source_ref:
            judgment_src.setdefault(src_key, []).append(j.source_ref)

    characters: list[Character] = []
    evidence: list[EvidenceLink] = []
    excluded_singletons = 0

    for (concept_id, cognate_class), members in sorted(class_members.items()):
        present_langs = set(members.keys())
        if drop_singletons and len(present_langs) < 2:
            excluded_singletons += 1
            continue

        states = [MISSING] * n_lang
        langs_with_concept = concept_langs.get(concept_id, set())
        for lid in languages:
            idx = lang_index[lid]
            if lid in present_langs:
                states[idx] = PRESENT
            elif lid in langs_with_concept:
                states[idx] = ABSENT
            else:
                states[idx] = MISSING

        loans = tuple(sorted(class_loans.get((concept_id, cognate_class), set())))
        characters.append(
            Character(
                concept_id=concept_id,
                cognate_class=cognate_class,
                states=states,
                is_partial=class_is_partial.get((concept_id, cognate_class), False),
                loan_languages=loans,
            )
        )

        # Evidence trail for present and absent cells (missing cells need none).
        for lid in languages:
            idx = lang_index[lid]
            st = states[idx]
            if st == MISSING:
                continue
            lex_ids = tuple(members.get(lid, []))
            srcs = tuple(judgment_src.get((concept_id, cognate_class, lid), []))
            evidence.append(
                EvidenceLink(
                    concept_id=concept_id,
                    cognate_class=cognate_class,
                    language_id=lid,
                    state=st,
                    lexeme_ids=lex_ids,
                    judgment_sources=srcs,
                )
            )

    matrix = CharacterMatrix(
        language_ids=languages,
        characters=characters,
        evidence=evidence,
        excluded_singletons=excluded_singletons,
    )
    return matrix
