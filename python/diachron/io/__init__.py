"""Data contracts: ingestion, validation, and the character transform.

Everything the models consume passes through this layer first. The guiding rule
is that the transformation from raw lexeme judgments to the character matrix is
explicit and reversible on paper: an audit sample can be traced from a cell in
the character matrix back to the lexemes and cognate judgments that produced it.
"""

from .types import (
    Concept,
    Language,
    Lexeme,
    CognateJudgment,
    Dataset,
    Attestation,
    DataError,
)

__all__ = [
    "Concept",
    "Language",
    "Lexeme",
    "CognateJudgment",
    "Dataset",
    "Attestation",
    "DataError",
]
