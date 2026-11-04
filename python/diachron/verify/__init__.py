"""Independent semantic verification of evidence bundles (gate G13).

The verifier is implemented separately from the summarizer that produces
evidence. It recomputes decisive quantities from raw records, performs bounded
replay using the exact saved engine, and rejects a bundle for a semantic reason
when a meaningful result or relationship has been altered. It never retrains,
regenerates, downloads a substitute, or trusts saved pass booleans.
"""

from .verifier import verify_evidence

__all__ = ["verify_evidence"]
