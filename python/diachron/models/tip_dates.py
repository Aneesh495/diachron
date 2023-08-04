"""Probability densities for dated terminal varieties."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .priors import lognormal_logpdf, normal_logpdf


@dataclass(frozen=True)
class TipDate:
    kind: str
    mean: float
    sd: float
    offset: float = 0.0

    def __post_init__(self) -> None:
        if (
            self.kind not in {"normal", "offset_lognormal"}
            or any(not math.isfinite(x) for x in (self.mean, self.sd, self.offset))
            or self.mean <= 0 or self.sd <= 0 or self.offset < 0
            or (self.kind == "normal" and self.offset != 0)
        ):
            raise ValueError("invalid dated-tip distribution")

    def center(self) -> float:
        if self.kind == "normal":
            return self.mean
        return self.offset + self.mean * math.exp(-self.sd ** 2 / 2)

    def log_density(self, age: float) -> float:
        if age <= 0 or not math.isfinite(age):
            return -math.inf
        if self.kind == "normal":
            # The calendar age is nonnegative; condition the normal accordingly.
            normalizer = 0.5 * math.erfc(-self.mean / (self.sd * math.sqrt(2.0)))
            return normal_logpdf(age, self.mean, self.sd) - math.log(normalizer)
        log_mean = math.log(self.mean) - self.sd ** 2 / 2
        return lognormal_logpdf(age - self.offset, log_mean, self.sd)

    def window(self) -> float:
        return self.sd if self.kind == "normal" else self.center() - self.offset
