"""Alternative-history registry (D24) and attestation inspector (D27)."""

from .scenarios import ScenarioRegistry, Scenario
from .inspector import AttestationInspector

__all__ = ["ScenarioRegistry", "Scenario", "AttestationInspector"]
