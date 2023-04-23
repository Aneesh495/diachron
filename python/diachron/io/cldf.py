"""CLDF ingestion into Diachron's typed tables.

CLDF (Cross-Linguistic Data Formats, https://cldf.clld.org/) packages lexical
data as CSV tables described by a JSON metadata file. A Wordlist module has a
LanguageTable, a ParameterTable (concepts), a FormTable (lexemes), and a
CognateTable (cognate judgments), each a CSV with CLDF-standard column
properties. This reader maps those tables onto `Concept`, `Language`, `Lexeme`,
and `CognateJudgment` records, carrying identifiers verbatim.

The reader is schema-driven: it reads the metadata JSON to find each table and
the column that plays each CLDF role (via the `propertyUrl`), rather than
assuming fixed column names, so it works across CLDF datasets that name their
columns differently. It records the dataset's declared license and, where
present, a DOI, because a shared schema does not grant one license to every
dataset.

Restricted datasets are never redistributed by Diachron. This reader ingests a
CLDF directory the user has acquired under its own license; the public release
ships the acquisition recipe, not the lexical records.
"""

from __future__ import annotations

import csv
import json
import os
from typing import Optional

from .types import (
    Concept,
    Language,
    Lexeme,
    CognateJudgment,
    Dataset,
    DataError,
)

# CLDF component property URLs that identify each table's role.
_WORDLIST_TABLES = {
    "LanguageTable": "http://cldf.clld.org/v1.0/terms.rdf#LanguageTable",
    "ParameterTable": "http://cldf.clld.org/v1.0/terms.rdf#ParameterTable",
    "FormTable": "http://cldf.clld.org/v1.0/terms.rdf#FormTable",
    "CognateTable": "http://cldf.clld.org/v1.0/terms.rdf#CognateTable",
}


def _find_column(table_schema: dict, property_suffix: str) -> Optional[str]:
    """Return the CSV column name bound to a CLDF property, or None."""
    for col in table_schema.get("tableSchema", {}).get("columns", []):
        prop = col.get("propertyUrl", "")
        if prop.endswith(property_suffix):
            return col["name"]
    return None


def _read_csv(base_dir: str, url: str) -> list[dict]:
    path = os.path.join(base_dir, url)
    if not os.path.exists(path):
        raise DataError(f"CLDF table file missing: {url}")
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_cldf_wordlist(metadata_path: str) -> Dataset:
    """Read a CLDF Wordlist directory into a Diachron `Dataset`.

    `metadata_path` points at the `*-metadata.json` file. The function locates
    the four Wordlist tables by their `dc:conformsTo` component URL and reads the
    role-bearing columns by property, so column naming differences do not break
    it.
    """
    with open(metadata_path, encoding="utf-8") as f:
        meta = json.load(f)
    base_dir = os.path.dirname(os.path.abspath(metadata_path))

    # Dataset-level license and identifier.
    license_ = meta.get("dc:license", "unknown")
    doi = None
    for ident in meta.get("dc:identifier", []) if isinstance(meta.get("dc:identifier"), list) else []:
        if isinstance(ident, str) and "doi.org" in ident:
            doi = ident
    name = meta.get("dc:title", os.path.basename(base_dir))
    version = str(meta.get("dc:hasVersion", meta.get("prov:wasGeneratedBy", "unknown")))

    tables: dict[str, dict] = {}
    for table in meta.get("tables", []):
        conforms = table.get("dc:conformsTo", "")
        for role, url in _WORDLIST_TABLES.items():
            if conforms == url:
                tables[role] = table

    if "FormTable" not in tables or "LanguageTable" not in tables:
        raise DataError("CLDF metadata missing a Form or Language table")

    # Languages.
    lang_table = tables["LanguageTable"]
    lang_id_col = _find_column(lang_table, "#id") or "ID"
    lang_name_col = _find_column(lang_table, "#name") or "Name"
    glotto_col = _find_column(lang_table, "#glottocode")
    lat_col = _find_column(lang_table, "#latitude")
    lon_col = _find_column(lang_table, "#longitude")
    languages: list[Language] = []
    for row in _read_csv(base_dir, lang_table["url"]):
        languages.append(
            Language(
                id=row[lang_id_col],
                name=row.get(lang_name_col, row[lang_id_col]),
                glottocode=row.get(glotto_col) if glotto_col else None,
                latitude=_float_or_none(row.get(lat_col)) if lat_col else None,
                longitude=_float_or_none(row.get(lon_col)) if lon_col else None,
            )
        )

    # Concepts (parameters). Some datasets omit a ParameterTable; derive concepts
    # from the forms in that case.
    concepts: list[Concept] = []
    concept_ids: set[str] = set()
    if "ParameterTable" in tables:
        par_table = tables["ParameterTable"]
        par_id = _find_column(par_table, "#id") or "ID"
        par_name = _find_column(par_table, "#name") or "Name"
        concepticon = _find_column(par_table, "#concepticonReference")
        for row in _read_csv(base_dir, par_table["url"]):
            concepts.append(
                Concept(
                    id=row[par_id],
                    gloss=row.get(par_name, row[par_id]),
                    concepticon_id=row.get(concepticon) if concepticon else None,
                )
            )
            concept_ids.add(row[par_id])

    # Forms (lexemes).
    form_table = tables["FormTable"]
    f_id = _find_column(form_table, "#id") or "ID"
    f_lang = _find_column(form_table, "#languageReference") or "Language_ID"
    f_param = _find_column(form_table, "#parameterReference") or "Parameter_ID"
    f_form = _find_column(form_table, "#form") or "Form"
    f_value = _find_column(form_table, "#value")
    f_source = _find_column(form_table, "#source")
    lexemes: list[Lexeme] = []
    for row in _read_csv(base_dir, form_table["url"]):
        cid = row[f_param]
        if cid not in concept_ids:
            # Derive a concept on first sight when no ParameterTable existed.
            concepts.append(Concept(id=cid, gloss=cid))
            concept_ids.add(cid)
        lexemes.append(
            Lexeme(
                id=row[f_id],
                language_id=row[f_lang],
                concept_id=cid,
                form=row.get(f_form, ""),
                original=row.get(f_value) if f_value else None,
                source_ref=row.get(f_source) if f_source else None,
            )
        )

    # Cognate judgments.
    judgments: list[CognateJudgment] = []
    if "CognateTable" in tables:
        cog_table = tables["CognateTable"]
        c_form = _find_column(cog_table, "#formReference") or "Form_ID"
        c_class = _find_column(cog_table, "#cognatesetReference") or "Cognateset_ID"
        c_source = _find_column(cog_table, "#source")
        form_concept = {lx.id: lx.concept_id for lx in lexemes}
        for row in _read_csv(base_dir, cog_table["url"]):
            form_id = row[c_form]
            concept_id = form_concept.get(form_id)
            if concept_id is None:
                continue
            judgments.append(
                CognateJudgment(
                    lexeme_id=form_id,
                    concept_id=concept_id,
                    cognate_class=row[c_class],
                    source_ref=row.get(c_source) if c_source else None,
                )
            )

    dataset = Dataset(
        name=name,
        version=version,
        license=license_,
        doi=doi,
        concepts=concepts,
        languages=languages,
        lexemes=lexemes,
        judgments=judgments,
        is_synthetic=False,
        source_url=metadata_path,
    )
    dataset.validate()
    return dataset


def _float_or_none(v) -> Optional[float]:
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
