"""Read a frozen, independently sourced EDH inscription cohort."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

from .freeze import verify_inputs
from .types import Attestation, DataError

_NS = {"t": "http://www.tei-c.org/ns/1.0"}
_XML_LANG = "{http://www.w3.org/XML/1998/namespace}lang"
_BP_EPOCH = 1950
_INSCRIPTION_LICENSE = "creativecommons.org/licenses/by-sa/4.0/"


def _catalogue_coordinates(path: Path) -> dict:
    """Index precise matches only, preserving ambiguity as a missing point."""
    points: dict[tuple[str, int], set[tuple[float, float]]] = defaultdict(set)
    features = json.loads(path.read_text(encoding="utf-8"))["features"]
    for feature in features:
        props = feature["properties"]
        name = props.get("ancient_findspot")
        uri = props.get("trismegistos_geo_uri")
        geometry = feature["geometry"]
        if not name or not uri or not geometry or geometry["type"] != "Point":
            continue
        tm = _trismegistos_id(uri)
        if tm is None:
            continue
        lon, lat = geometry["coordinates"]
        if -90 <= lat <= 90 and -180 <= lon <= 180:
            points[(name.casefold(), tm)].add((lat, lon))
    return {key: next(iter(coords)) for key, coords in points.items() if len(coords) == 1}


def _trismegistos_id(uri: str) -> int | None:
    match = re.search(r"trismegistos\.org/place/(\d+)$", uri)
    return int(match[1]) if match else None


def _bp(year: int) -> int:
    # EDH signed dates use historical years without a year zero.
    return _BP_EPOCH - (year + (year < 0))


def _read_inscription(path: Path, places: dict, latin_tip: str) -> tuple[list[Attestation], str]:
    raw = path.read_bytes()
    if b"<!DOCTYPE" in raw.upper():
        raise DataError(f"unexpected document type in {path.name}")
    root = ET.fromstring(raw)
    identifier = root.findtext(".//t:idno[@type='localID']", namespaces=_NS)
    if identifier != path.stem or not re.fullmatch(r"HD\d{6}", identifier):
        raise DataError(f"EDH identifier disagrees with filename: {path.name}")
    licence = root.find(".//t:availability/t:licence", _NS)
    if licence is None or _INSCRIPTION_LICENSE not in licence.get("target", ""):
        raise DataError(f"EDH licence missing in {path.name}")
    languages = {
        edition.get(_XML_LANG) for edition in root.findall(".//t:div[@type='edition']", _NS)
    }
    if languages != {"la"}:
        return [], "language"
    bibliography = tuple(
        " ".join("".join(ref.itertext()).split())
        for ref in root.findall(".//t:div[@type='bibliography']//t:bibl", _NS)
        if "".join(ref.itertext()).strip()
    )
    if not bibliography:
        return [], "bibliography"
    origin = root.find(".//t:history/t:origin", _NS)
    if origin is None:
        return [], "origin"
    localities = [
        name for name in origin.findall("t:origPlace/t:placeName", _NS)
        if name.text and name.text.strip() and not name.get("type")
    ]
    if len(localities) != 1:
        return [], "locality"
    locality = localities[0]
    region = locality.text.strip()
    tm = _trismegistos_id(locality.get("ref", ""))
    point = places.get((region.casefold(), tm)) if tm is not None else None
    latitude, longitude = point or (None, None)

    intervals = []
    for date in origin.findall("t:origDate", _NS):
        low = date.get("notBefore-custom")
        high = date.get("notAfter-custom")
        if not low or not high:
            continue
        start = re.fullmatch(r"(-?\d{4})(?:-\d{2}(?:-\d{2})?)?", low)
        end = re.fullmatch(r"(-?\d{4})(?:-\d{2}(?:-\d{2})?)?", high)
        if start is None or end is None:
            raise DataError(f"unrecognized EpiDoc date range in {path.name}: {low}, {high}")
        lo, hi = int(start[1]), int(end[1])
        if lo == 0 or hi == 0 or hi < lo or hi > _BP_EPOCH:
            continue
        interval = (_bp(hi), _bp(lo), low, high)
        if interval not in intervals:
            intervals.append(interval)
    if not intervals:
        return [], "date"

    ids = [
        f"edh:{identifier}" + (f":date-{i + 1}" if len(intervals) > 1 else "")
        for i in range(len(intervals))
    ]
    return [
        Attestation(
            id=obs_id, language_id=latin_tip, age_low=lo, age_high=hi,
            source_ref=f"https://edh.ub.uni-heidelberg.de/edh/inschrift/{identifier}",
            page="EpiDoc teiHeader/fileDesc/sourceDesc/msDesc/history/origin",
            region=region, latitude=latitude, longitude=longitude,
            language_assignment_confidence=None, observation_type="inscription",
            conflicting_with=tuple(alt for alt in ids if alt != obs_id),
            bibliography=bibliography, location_kind="reported_origin",
            source_date_start=source_start, source_date_end=source_end,
            location_uncertainty_note=(
                "EDH locality reference point. Object origin precision is unquantified."
                if point else "Source locality only. No unique matching gazetteer point."
            ),
            language_assignment_note=(
                "EDH edition tagged Latin; mapped to the IE-CoR Latin variety for indexing."
            ),
        )
        for obs_id, (lo, hi, source_start, source_end) in zip(ids, intervals)
    ], "included"


def load_edh(source: Path, geography: Path) -> dict:
    """Return source records and exclusions from the frozen selection."""
    manifest = json.loads(Path("data/manifests/edh.json").read_text(encoding="utf-8"))
    verify_inputs(source, manifest["files"])
    verify_inputs(geography.parent, {geography.name: manifest["geography"]})
    places = _catalogue_coordinates(geography)
    attestations = []
    dispositions: Counter[str] = Counter()
    for name in manifest["files"]:
        observations, reason = _read_inscription(source / name, places, manifest["iecor_latin_id"])
        attestations.extend(observations)
        dispositions[reason] += 1
    if dispositions["included"] < 100:
        raise DataError("EDH cohort has fewer than 100 independently identified inscriptions")
    return {
        "schema_version": 1,
        "source": manifest["url"],
        "source_commit": manifest["git_commit"],
        "licence": manifest["license"],
        "bp_epoch": _BP_EPOCH,
        "is_synthetic": False,
        "dispositions": dict(dispositions),
        "n_attestations": len(attestations),
        "n_with_catalogue_point": sum(a.latitude is not None for a in attestations),
        "attestations": [asdict(a) for a in attestations],
    }
