"""Geographic models: region CTMC and spherical diffusion.

Two geographic observation models are offered and compared rather than assumed.
The region model treats location as a discrete latent state over coarse named
regions and reuses the discrete pruning likelihood; the spherical model treats
location as a continuous diffusing coordinate. Both distinguish proximity,
documented contact, and hypothetical contact, and neither applies a modern
coastline or political border to all historical periods.
"""

from .regions import RegionModel, assign_regions, region_log_likelihood
from .contact import ContactKind, ContactEdge, classify_contact

__all__ = [
    "RegionModel",
    "assign_regions",
    "region_log_likelihood",
    "ContactKind",
    "ContactEdge",
    "classify_contact",
]
