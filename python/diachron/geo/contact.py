"""Contact classification for the borrowing/geography interface.

The brief insists that proximity, documented contact, and hypothetical contact be
kept distinct. A candidate borrowing edge between two lineages is labeled by the
kind of evidence that supports it, and the label travels with the edge into the
borrowing model and the UI so a merely-nearby pair is never presented as a
documented contact.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from diachron import _core


class ContactKind(str, Enum):
    PROXIMITY = "proximity"            # near in space, no documentation
    DOCUMENTED = "documented"          # attested loan or recorded contact
    HYPOTHETICAL = "hypothetical"      # posited by a hypothesis, not documented


@dataclass
class ContactEdge:
    donor_taxon: int
    recipient_taxon: int
    kind: ContactKind
    distance_km: float
    note: str = ""


def classify_contact(
    donor_taxon: int,
    recipient_taxon: int,
    donor_coord: tuple[float, float],
    recipient_coord: tuple[float, float],
    documented: bool,
    hypothesized: bool,
    proximity_km: float = 600.0,
) -> ContactEdge:
    """Classify a candidate contact edge by the strongest evidence it carries.

    Documentation outranks a hypothesis, which outranks mere proximity. An edge
    that is neither documented, hypothesized, nor within the proximity threshold
    is still returned as proximity-kind with its true distance, so the caller can
    filter it out; the classification never upgrades evidence it does not have.
    """
    d = _core.haversine_km(
        donor_coord[0], donor_coord[1], recipient_coord[0], recipient_coord[1]
    )
    if documented:
        kind = ContactKind.DOCUMENTED
    elif hypothesized:
        kind = ContactKind.HYPOTHETICAL
    else:
        kind = ContactKind.PROXIMITY
    note = ""
    if kind is ContactKind.PROXIMITY and d > proximity_km:
        note = f"beyond proximity threshold ({proximity_km} km); weak candidate"
    return ContactEdge(
        donor_taxon=donor_taxon,
        recipient_taxon=recipient_taxon,
        kind=kind,
        distance_km=d,
        note=note,
    )
