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
from pathlib import Path
from unicodedata import normalize
from typing import Optional

from .types import (
    Concept,
    Language,
    Lexeme,
    CognateJudgment,
    CognateClassRecord,
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


def _read_csv(base_dir: str, url: str):
    base = Path(base_dir).resolve()
    path = (base / url).resolve()
    if not path.is_relative_to(base) or not path.is_file():
        raise DataError(f"CLDF table must be a local file inside the dataset: {url}")
    if path.stat().st_size > 128 * 1024 * 1024:
        raise DataError(f"CLDF table exceeds 128 MiB: {url}")
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = csv.DictReader(f)
        if not rows.fieldnames or len(set(rows.fieldnames)) != len(rows.fieldnames):
            raise DataError(f"CLDF table has an invalid header: {url}")
        for row in rows:
            if None in row or None in row.values():
                raise DataError(f"CLDF table has a malformed row at {url}:{rows.line_num}")
            yield row


def load_cldf_wordlist(
    metadata_path: str, *, profile: str | None = None,
) -> Dataset:
    """Read a CLDF Wordlist directory into a Diachron `Dataset`.

    `metadata_path` points at the `*-metadata.json` file. The function locates
    the four Wordlist tables by their `dc:conformsTo` component URL and reads the
    role-bearing columns by property, so column naming differences do not break
    it.
    """
    if profile not in (None, "iecor", "abvd"):
        raise DataError(f"unknown CLDF source profile: {profile}")
    if Path(metadata_path).stat().st_size > 4 * 1024 * 1024:
        raise DataError("CLDF metadata exceeds 4 MiB")
    with open(metadata_path, encoding="utf-8") as f:
        meta = json.load(f)
    base_dir = os.path.dirname(os.path.abspath(metadata_path))

    # Dataset-level license and identifier.
    license_ = meta.get("dc:license", "unknown")
    identifiers = meta.get("dc:identifier", [])
    if isinstance(identifiers, str):
        identifiers = [identifiers]
    doi = next(
        (ident for ident in identifiers if isinstance(ident, str) and "doi.org/" in ident),
        None,
    )
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
        date = row.get("distribution") if profile == "iecor" else None
        parameters = {}
        if date == "Normal":
            parameters = dict(
                date_distribution="normal", date_mean=float(row["normalMean"]),
                date_sd=float(row["normalStDev"]), date_offset=0.0,
            )
        elif date == "Offset log normal":
            parameters = dict(
                date_distribution="offset_lognormal", date_mean=float(row["logNormalMean"]),
                date_sd=float(row["logNormalStDev"]), date_offset=float(row["logNormalOffset"]),
            )
        elif date:
            raise DataError(f"unknown source date distribution: {date}")
        languages.append(
            Language(
                id=row[lang_id_col],
                name=row.get(lang_name_col, row[lang_id_col]),
                glottocode=row.get(glotto_col) if glotto_col else None,
                latitude=_float_or_none(row.get(lat_col)) if lat_col else None,
                longitude=_float_or_none(row.get(lon_col)) if lon_col else None,
                family=(row.get("Family") or None) if profile else None,
                is_historical=_flag(row.get("historical")) is True if profile == "iecor" else False,
                source_classification=(row.get("Clade") or None) if profile == "iecor" else None,
                **parameters,
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
            if "ParameterTable" in tables:
                raise DataError(f"form {row[f_id]} references unknown concept {cid}")
            concepts.append(Concept(id=cid, gloss=cid))
            concept_ids.add(cid)
        raw_value = (row.get(f_value) or None) if f_value else None
        native = (row.get("native_script") or None) if profile == "iecor" else None
        lexemes.append(
            Lexeme(
                id=row[f_id],
                language_id=row[f_lang],
                concept_id=cid,
                form=normalize("NFC", row.get(f_form, "")),
                original=native or raw_value,
                source_value=raw_value,
                transliteration_source=f"{form_table['url']}:{f_form}" if native else None,
                loan_annotation=_flag(row.get("Loan")) if profile else None,
                source_ref=(row.get(f_source) or None) if f_source else None,
            )
        )

    # Cognate judgments.
    judgments: list[CognateJudgment] = []
    if "CognateTable" in tables:
        cog_table = tables["CognateTable"]
        c_form = _find_column(cog_table, "#formReference") or "Form_ID"
        c_class = _find_column(cog_table, "#cognatesetReference") or "Cognateset_ID"
        c_source = _find_column(cog_table, "#source")
        c_slice = _find_column(cog_table, "#segmentSlice")
        c_id = _find_column(cog_table, "#id") or "ID"
        forms = {lx.id: lx for lx in lexemes}
        loans = {}
        if profile == "iecor":
            loans = {row["Cognateset_ID"]: row for row in _read_csv(base_dir, "loans.csv")}
        judgment_ids = set()
        for row in _read_csv(base_dir, cog_table["url"]):
            form_id = row[c_form]
            form = forms.get(form_id)
            if form is None:
                raise DataError(f"cognate judgment references unknown form {form_id}")
            identifier = row.get(c_id)
            if not identifier or identifier in judgment_ids:
                raise DataError(f"missing or duplicate cognate judgment id: {identifier}")
            judgment_ids.add(identifier)
            cognate_class = row[c_class]
            segment_slice = (row.get(c_slice) or None) if c_slice else None
            loan = loans.get(cognate_class)
            judgments.append(
                CognateJudgment(
                    lexeme_id=form_id,
                    concept_id=form.concept_id,
                    cognate_class=cognate_class,
                    source_ref=(row.get(c_source) or None) if c_source else None,
                    id=identifier,
                    doubtful=_flag(row.get("Doubt")) if profile else None,
                    detection_method=(row.get("Cognate_Detection_Method") or None) if profile else None,
                    segment_slice=segment_slice,
                    is_partial=segment_slice is not None,
                    is_loan=loan is not None or form.loan_annotation is True,
                    loan_event=cognate_class if loan is not None else None,
                    loan_source=(loan.get("SourceCognateset_ID") or None) if loan else None,
                )
            )

    classes = []
    if profile == "iecor":
        for row in _read_csv(base_dir, "cognatesets.csv"):
            scale = row.get("proposedAsCognateToScale")
            classes.append(CognateClassRecord(
                id=row["ID"], root_form=row.get("Root_Form") or None,
                root_language=row.get("Root_Language") or None,
                source_ref=row.get("Source") or None,
                proposed_cognate_id=row.get("proposedAsCognateTo_pk") or None,
                proposed_scale=int(scale) if scale else None,
                superset_id=row.get("supersetid") or None,
                parallel_derivation=_flag(row.get("parallelDerivation")),
            ))
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
        cognate_classes=classes,
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


def _flag(value: str | None) -> bool | None:
    if not value:
        return None
    if value.lower() not in ("true", "false", "1", "0"):
        raise DataError(f"invalid source annotation flag: {value!r}")
    return value.lower() in ("true", "1")
