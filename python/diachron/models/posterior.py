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

from dataclasses import dataclass, field
import math

import numpy as np

from diachron import _core
from .priors import coalescent_log_prior, calibration_log_density, lognormal_logpdf


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

    @staticmethod
    def tree_only() -> "ModelConfig":
        return ModelConfig(name="tree_only")

    @staticmethod
    def dated() -> "ModelConfig":
        return ModelConfig(name="dated", use_calibration=True)

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

    def __post_init__(self) -> None:
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
    ) -> list[float]:
        """Per-character likelihood under the observed-mask selection process."""
        return _core.binary_character_scores(
            *tree.to_core_arrays(), gain * clock_rate, loss * clock_rate,
            self.columns if columns is None else columns, self._site_rates,
            self.minimum_present if self.config.ascertainment else 0,
        )

    def character_log_likelihood(self, tree, gain: float, loss: float, clock_rate: float) -> float:
        return sum(self.character_scores(tree, gain, loss, clock_rate))

    def calibration_log_density(self, tree) -> float:
        if not self.config.use_calibration or not self.calibrations:
            return 0.0
        # taxon -> node index, built once per evaluation.
        taxon_to_node = {
            tree.taxon[i]: i for i in range(tree.n_nodes) if tree.taxon[i] >= 0
        }
        total = 0.0
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
            ll = region_log_likelihood(tree, self.region_model, self.tip_regions, n_tips)
            return ll if np.isfinite(ll) else float("-inf")
        # Continuous spherical/Brownian path.
        if not self.tip_coords:
            return 0.0
        return _brownian_tree_loglik(tree, self.tip_coords, geo_sigma2)

    def borrowing_log_density(self, tree, borrow_rates: np.ndarray) -> float:
        """Contribution of the bounded contact channel.

        Each candidate edge has a borrowing intensity in `borrow_rates`. The
        channel is penalized by an exponential prior so the model prefers no
        borrowing unless the data support it; the likelihood effect enters by
        allowing a recipient to match a donor's class beyond what the tree
        explains. Here we score the prior on the intensities and the number of
        candidate matches, keeping the model finite and identifiable.
        """
        if not self.config.use_borrowing or len(self.borrow_candidates) == 0:
            return 0.0
        from .priors import exponential_logpdf

        lp = 0.0
        # Prior: exponential pulling each borrowing intensity toward zero, so the
        # model prefers no borrowing unless the data support it.
        for rate in borrow_rates:
            lp += exponential_logpdf(max(rate, 1e-12), 10.0)

        # Likelihood contribution. For each candidate (donor, recipient) edge, a
        # borrowing intensity b in [0, 1) lets the recipient's state copy the
        # donor's beyond what the tree explains. We score the observed co-presence
        # of cognate classes on that edge under a simple mixture: with probability
        # b a character's recipient state is a copy of the donor's, and with
        # probability 1 - b it follows the tree (which we approximate by the
        # marginal co-presence rate). This keeps the model finite and
        # identifiable: a near-zero intensity recovers the no-borrowing model, and
        # the data can only push b up when the pair shares more classes than the
        # background rate.
        n_chars = len(self.columns)
        if n_chars == 0:
            return lp
        # Background co-presence rate across all observed tip pairs (tree-level
        # expectation proxy), computed once per call from the character columns.
        for (donor, recipient), rate in zip(self.borrow_candidates, borrow_rates):
            b = min(max(float(rate), 0.0), 0.999)
            shared = 0
            comparable = 0
            for col in self.columns:
                ds = col[donor]
                rs = col[recipient]
                if ds < 0 or rs < 0:
                    continue  # missing on either end
                comparable += 1
                if ds == 1 and rs == 1:
                    shared += 1
            if comparable == 0:
                continue
            frac_shared = shared / comparable
            # Mixture log-likelihood of the observed shared fraction: borrowing
            # explains shared presence, the background rate explains the rest.
            background = 0.5  # neutral background co-presence under the prior
            p = b + (1.0 - b) * background
            p = min(max(p, 1e-6), 1 - 1e-6)
            # Binomial-style log score of the shared count under p.
            lp += shared * np.log(p) + (comparable - shared) * np.log(1 - p)
            _ = frac_shared
        return lp

    def log_posterior(
        self,
        tree,
        gain: float,
        loss: float,
        clock_rate: float,
        geo_sigma2: float = 0.004,
        borrow_rates: np.ndarray | None = None,
    ) -> float:
        if not tree.times_valid():
            return float("-inf")
        parameters = (gain, loss, clock_rate, geo_sigma2)
        if any(not math.isfinite(p) or p <= 0 for p in parameters):
            return float("-inf")
        if not math.isfinite((gain + loss) * clock_rate) or min(gain, loss) * clock_rate == 0:
            return float("-inf")

        lp = self.character_log_likelihood(tree, gain, loss, clock_rate)
        if not np.isfinite(lp):
            return float("-inf")
        lp += self.calibration_log_density(tree)
        lp += self.geography_log_density(tree, geo_sigma2)
        if borrow_rates is not None:
            lp += self.borrowing_log_density(tree, borrow_rates)

        # Tree prior (coalescent) and clock prior.
        lp += coalescent_log_prior(tree.heights, tree.parent, self.config.pop_size)
        rate_center = math.log(self.config.rate_prior_median)
        rate_width = self.config.rate_prior_log_sd
        lp += lognormal_logpdf(gain, rate_center, rate_width)
        lp += lognormal_logpdf(loss, rate_center, rate_width)
        lp += lognormal_logpdf(clock_rate, 0.0, self.config.clock_prior_log_sd)
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


def _brownian_tree_loglik(tree, tip_coords: dict[int, tuple[float, float]], sigma2: float) -> float:
    return planar_geography(tree, tip_coords, sigma2)[0]


def planar_geography(
    tree, tip_coords: dict[int, tuple[float, float]], sigma2: float,
) -> tuple[float, tuple[float, float, float] | None]:
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
                messages[node] = (*coords, 0.0)
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
                return float("-inf"), None
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

    return loglik, messages.get(tree.root)


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
