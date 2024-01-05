"""The composed log posterior and the model configurations compared in the study.

`ModelConfig` names which terms are active, giving the four required models:

* **tree_only**: character likelihood + tree prior + strict clock prior. No
  calibration, so only relative times are identified (the root age floats).
* **dated**: adds node-age calibrations, making the time scale identifiable.
* **geography**: adds a geographic likelihood over tip locations and latent
  internal locations under a diffusion model.
* **borrowing**: adds a bounded contact channel that lets a small set of
  candidate edges donate cognate classes, so tree-versus-borrowing conflict can
  be examined.

The observation models are directly comparable: every model scores the same
character matrix with the same ascertainment correction, so differences in fit
come from the added structure, not from changing what counts as data.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

import numpy as np

from diachron import _core
from .tip_dates import TipDate
from .constraints import validate_clades, satisfies_clades
from .priors import coalescent_log_prior, calibration_log_density, exponential_logpdf, lognormal_logpdf


@dataclass
class Calibration:
    """A node-age calibration tying an internal node to a dated interval."""

    # Tip taxa whose most recent common ancestor the calibration constrains.
    clade_taxa: tuple[int, ...]
    age_low: float
    age_high: float
    kind: str = "uniform"
    source_ref: str = ""


@dataclass
class ModelConfig:
    name: str
    use_calibration: bool = False
    use_geography: bool = False
    use_borrowing: bool = False
    pop_size: float = 3000.0  # coalescent effective scale in years
    ascertainment: bool = True
    # Among-character rate heterogeneity. None disables it (single rate); a
    # positive alpha enables a discrete-gamma mixture with that shape and
    # `rate_categories` categories.
    rate_alpha: float | None = None
    rate_categories: int = 4
    rate_prior_median: float = 0.001  # events per year before the clock multiplier
    rate_prior_log_sd: float = 1.5
    clock_prior_log_sd: float = 0.5  # dimensionless multiplier centered at one
    diffusion_prior_median: float = 0.004  # degree squared / year
    diffusion_prior_log_sd: float = 2.0
    borrowing_prior_rate: float = 10.0
    relaxed_clock_sigma: float = 0.0
    geo_region_count: int = 0
    geography_error_degrees: float = 0.0

    @staticmethod
    def tree_only() -> "ModelConfig":
        return ModelConfig(name="tree_only")

    @staticmethod
    def dated() -> "ModelConfig":
        return ModelConfig(name="dated", use_calibration=True)

    @staticmethod
    def relaxed() -> "ModelConfig":
        return ModelConfig(name="relaxed", use_calibration=True, relaxed_clock_sigma=0.25)

    @staticmethod
    def geography() -> "ModelConfig":
        return ModelConfig(name="geography", use_calibration=True, use_geography=True)

    @staticmethod
    def borrowing() -> "ModelConfig":
        return ModelConfig(
            name="borrowing", use_calibration=True, use_borrowing=True
        )


@dataclass
class PosteriorModel:
    """Everything needed to evaluate the log posterior of a state.

    The data (character columns, tip coordinates, calibrations) are fixed; the
    sampler mutates the tree and the scalar parameters and asks this object for
    the resulting unnormalized log posterior.
    """

    config: ModelConfig
    columns: list[list[int]]           # variable character columns, binding convention
    n_states: int = 2
    calibrations: list[Calibration] = field(default_factory=list)
    tip_coords: dict[int, tuple[float, float]] = field(default_factory=dict)
    # Bounded candidate borrowing edges as (donor_taxon, recipient_taxon).
    borrow_candidates: list[tuple[int, int]] = field(default_factory=list)
    # Optional region model; when set, geography uses the discrete region CTMC
    # instead of the continuous Brownian likelihood.
    region_model: object = None
    tip_regions: dict[int, int] = field(default_factory=dict)
    n_tips: int = 0

    # Two when singleton positives were filtered, one for variable-only data.
    minimum_present: int = 1
    tip_dates: dict[int, TipDate] = field(default_factory=dict)
    required_clades: list[list[int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._clade_masks = validate_clades(self.required_clades, self.n_tips)
        if not math.isfinite(self.config.relaxed_clock_sigma) or self.config.relaxed_clock_sigma < 0.0:
            raise ValueError("clock dispersion must be finite and nonnegative")
        error = self.config.geography_error_degrees
        if not math.isfinite(error) or error < 0:
            raise ValueError("geographic observation error must be finite and nonnegative")
        alpha = self.config.rate_alpha
        if alpha is not None and (
            not np.isfinite(alpha) or alpha <= 0 or not 1 <= self.config.rate_categories <= 32
        ):
            raise ValueError("gamma shape must be positive with 1 to 32 rate categories")
        self._site_rates = (
            [1.0] if alpha is None
            else _core.discrete_gamma_rates(alpha, self.config.rate_categories)
        )

    def character_scores(
        self, tree, gain: float, loss: float, clock_rate: float,
        columns: list[list[int]] | None = None,
        borrow_rates: list[float] | None = None,
        clock_multipliers: list[float] | None = None,
    ) -> list[float]:
        """Per-character likelihood under the observed-mask selection process."""
        contacts = []
        if self.config.use_borrowing and borrow_rates is not None:
            if len(borrow_rates) != len(self.borrow_candidates):
                raise ValueError("borrowing intensity count must match the fixed candidate network")
            contacts = [
                (donor, recipient, intensity)
                for (donor, recipient), intensity in zip(self.borrow_candidates, borrow_rates)
                if intensity > 0.0
            ]
        parent, lengths, heights, taxa, root = tree.to_core_arrays()
        if clock_multipliers:
            lengths = [t * multiplier for t, multiplier in zip(lengths, clock_multipliers)]
        return _core.binary_character_scores(
            parent, lengths, heights, taxa, root, gain * clock_rate, loss * clock_rate,
            self.columns if columns is None else columns, self._site_rates,
            self.minimum_present if self.config.ascertainment else 0, contacts,
        )

    def character_log_likelihood(self, tree, gain: float, loss: float, clock_rate: float) -> float:
        return sum(self.character_scores(tree, gain, loss, clock_rate))

    def calibration_log_density(self, tree) -> float:
        if not self.config.use_calibration:
            return 0.0
        # taxon -> node index, built once per evaluation.
        taxon_to_node = {
            tree.taxon[i]: i for i in range(tree.n_nodes) if tree.taxon[i] >= 0
        }
        total = sum(
            date.log_density(tree.heights[taxon_to_node[taxon]])
            for taxon, date in self.tip_dates.items()
        )
        need_mrca = False
        for cal in self.calibrations:
            if len(cal.clade_taxa) == 1:
                # Single-taxon calibration from a dated attestation: the
                # language's lineage existed at the attested age, so the node
                # where that lineage branches from its sibling (the tip's parent)
                # must be at least that old. We score the parent node's age
                # against the interval, which is what carries time-scale
                # information; a contemporary tip at age 0 would carry none.
                tip = taxon_to_node.get(cal.clade_taxa[0])
                if tip is None:
                    return float("-inf")
                parent = tree.parent[tip]
                node = parent if parent >= 0 else tip
                total += calibration_log_density(
                    tree.heights[node], cal.age_low, cal.age_high, cal.kind
                )
            else:
                need_mrca = True
        if need_mrca:
            n_tips = len(taxon_to_node)
            for cal in self.calibrations:
                if len(cal.clade_taxa) <= 1:
                    continue
                if len(cal.clade_taxa) == n_tips:
                    # MRCA of every tip is the root; skip the ancestor walk.
                    node = tree.root
                else:
                    node = _mrca(tree, cal.clade_taxa, taxon_to_node)
                if node is None:
                    return float("-inf")
                total += calibration_log_density(
                    tree.heights[node], cal.age_low, cal.age_high, cal.kind
                )
        return total

    def geography_log_density(self, tree, geo_sigma2: float) -> float:
        """Geographic term under planar contrasts or an attached region model.

        The planar model uses raw latitude and longitude as a local approximation.
        Its degree-squared diffusion scale does not describe global movement.
        """
        if not self.config.use_geography:
            return 0.0
        # Region CTMC path: a discrete latent-region likelihood over the tree,
        # used when a region model is attached. This is the coarse-evidence model.
        if self.region_model is not None and self.tip_regions:
            from ..geo.regions import region_log_likelihood

            n_tips = self.n_tips or (max(self.tip_regions) + 1)
            region = replace(
                self.region_model, base_rate=self.region_model.base_rate * geo_sigma2 / 0.004,
            )
            ll = region_log_likelihood(tree, region, self.tip_regions, n_tips)
            return ll if np.isfinite(ll) else float("-inf")
        # Continuous spherical/Brownian path.
        if not self.tip_coords:
            return 0.0
        return _brownian_tree_loglik(
            tree, self.tip_coords, geo_sigma2, self.config.geography_error_degrees,
        )

    def log_posterior(
        self,
        tree,
        gain: float,
        loss: float,
        clock_rate: float,
        geo_sigma2: float = 0.004,
        borrow_rates: list[float] | None = None,
        clock_multipliers: list[float] | None = None,
    ) -> float:
        if not tree.times_valid() or not satisfies_clades(tree, self._clade_masks):
            return float("-inf")
        parameters = (gain, loss, clock_rate, geo_sigma2)
        if any(not math.isfinite(p) or p <= 0 for p in parameters):
            return float("-inf")
        if not math.isfinite((gain + loss) * clock_rate) or min(gain, loss) * clock_rate == 0:
            return float("-inf")

        if self.config.use_borrowing:
            if borrow_rates is None or len(borrow_rates) != len(self.borrow_candidates):
                return float("-inf")
            if any(not math.isfinite(b) or b < 0 for b in borrow_rates):
                return float("-inf")
        sigma = self.config.relaxed_clock_sigma
        if sigma > 0:
            if (
                clock_multipliers is None or len(clock_multipliers) != tree.n_nodes
                or any(not math.isfinite(m) or m <= 0 for m in clock_multipliers)
                or clock_multipliers[tree.root] != 1.0
            ):
                return float("-inf")
        elif clock_multipliers:
            return float("-inf")
        lp = sum(self.character_scores(
            tree, gain, loss, clock_rate, borrow_rates=borrow_rates,
            clock_multipliers=clock_multipliers,
        ))
        if not np.isfinite(lp):
            return float("-inf")
        lp += self.calibration_log_density(tree)
        lp += self.geography_log_density(tree, geo_sigma2)
        if self.config.use_borrowing:
            lp += sum(exponential_logpdf(rate, self.config.borrowing_prior_rate) for rate in borrow_rates)

        # Tree prior (coalescent) and clock prior.
        lp += coalescent_log_prior(tree.heights, tree.parent, self.config.pop_size)
        rate_center = math.log(self.config.rate_prior_median)
        rate_width = self.config.rate_prior_log_sd
        lp += lognormal_logpdf(gain, rate_center, rate_width)
        lp += lognormal_logpdf(loss, rate_center, rate_width)
        lp += lognormal_logpdf(clock_rate, 0.0, self.config.clock_prior_log_sd)
        if sigma > 0:
            lp += sum(
                lognormal_logpdf(multiplier, -sigma * sigma / 2, sigma)
                for node, multiplier in enumerate(clock_multipliers) if node != tree.root
            )
        if self.config.use_geography:
            lp += lognormal_logpdf(
                geo_sigma2, math.log(self.config.diffusion_prior_median),
                self.config.diffusion_prior_log_sd,
            )
        return lp


def _mrca(tree, taxa: tuple[int, ...], taxon_to_node: dict[int, int]) -> int | None:
    """Most recent common ancestor node of a set of tip taxa."""
    try:
        nodes = [taxon_to_node[t] for t in taxa]
    except KeyError:
        return None

    def ancestors(node: int) -> list[int]:
        path = [node]
        while tree.parent[node] >= 0:
            node = tree.parent[node]
            path.append(node)
        return path

    common = set(ancestors(nodes[0]))
    for n in nodes[1:]:
        common &= set(ancestors(n))
    if not common:
        return None
    # The MRCA is the youngest common ancestor.
    return min(common, key=lambda n: tree.heights[n])


def _brownian_tree_loglik(
    tree, tip_coords: dict[int, tuple[float, float]], sigma2: float, error: float = 0.0,
) -> float:
    return _planar_messages(tree, tip_coords, sigma2, error)[0]


def _planar_messages(
    tree, tip_coords: dict[int, tuple[float, float]], sigma2: float, error: float,
) -> tuple[float, dict[int, tuple[float, float, float]]]:
    """Planar Gaussian contrasts over observed latitude and longitude pairs.

    A message is a conditional mean and its variance. Missing subtrees
    contribute no message. Combining children sequentially is equivalent to
    integrating a Gaussian common ancestor and supports multifurcations.
    The likelihood is defined on contrasts, with a translation-invariant root.
    """
    children = tree.children()
    durations = tree.branch_lengths()
    messages = {}
    loglik = 0.0
    for node in _postorder(tree):
        if tree.is_tip(node):
            coords = tip_coords.get(tree.taxon[node])
            if coords is not None:
                messages[node] = (*coords, error ** 2)
            continue

        combined = None
        for child in children[node]:
            if child not in messages:
                continue
            lat, lon, variance = messages[child]
            variance += sigma2 * durations[child]
            if combined is None:
                combined = (lat, lon, variance)
                continue

            prev_lat, prev_lon, prev_var = combined
            contrast_var = prev_var + variance
            if contrast_var <= 0.0:
                return float("-inf"), {}
            distance2 = (lat - prev_lat) ** 2 + (lon - prev_lon) ** 2
            loglik -= 0.5 * distance2 / contrast_var + math.log(2 * math.pi * contrast_var)
            weight = prev_var / contrast_var
            combined = (
                prev_lat + weight * (lat - prev_lat),
                prev_lon + weight * (lon - prev_lon),
                weight * variance,
            )
        if combined is not None:
            messages[node] = combined

    return loglik, messages


def planar_location_marginals(
    tree, tip_coords: dict[int, tuple[float, float]], sigma2: float, error: float = 0.0,
) -> dict[int, tuple[float, float, float]]:
    """Condition every node's location on the complete set of scored tips."""
    _, upward = _planar_messages(tree, tip_coords, sigma2, error)
    if tree.root not in upward:
        return {}
    posterior = {tree.root: upward[tree.root]}
    durations = tree.branch_lengths()
    children = tree.children()
    order = [tree.root]
    for parent in order:
        order.extend(children[parent])
        parent_lat, parent_lon, parent_var = posterior[parent]
        for node in children[parent]:
            edge_var = sigma2 * durations[node]
            if node not in upward:
                posterior[node] = (parent_lat, parent_lon, parent_var + edge_var)
                continue
            local_lat, local_lon, local_var = upward[node]
            combined = local_var + edge_var
            if combined == 0:
                posterior[node] = upward[node]
                continue
            parent_weight = local_var / combined
            local_weight = edge_var / combined
            posterior[node] = (
                parent_weight * parent_lat + local_weight * local_lat,
                parent_weight * parent_lon + local_weight * local_lon,
                parent_weight * (edge_var + parent_weight * parent_var),
            )
    return posterior


def _postorder(tree) -> list[int]:
    ch = tree.children()
    order: list[int] = []
    stack = [(tree.root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            order.append(node)
        else:
            stack.append((node, True))
            for c in ch[node]:
                stack.append((c, False))
    return order
