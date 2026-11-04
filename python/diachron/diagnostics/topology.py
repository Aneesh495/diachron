"""Topology mixing diagnostics.

Scalar ESS does not tell us whether the discrete tree topology mixed, so we
summarize topology sampling directly: how often each distinct topology appears,
how often the sampled topology changes between consecutive draws (a proxy for
mixing speed), and whether independent chains agree on the posterior over clades.

Topologies are identified by their set of clades (bipartitions induced by each
internal node), which is invariant to node-index relabeling, so two trees with
the same shape and tip labels are recognized as the same topology.
"""

from __future__ import annotations

from collections import Counter
import re


def _clade_signature(newick: str) -> frozenset[frozenset[str]]:
    """Extract the set of tip-label clades from a Newick string.

    We parse the nested parenthesization and record, for each internal node, the
    set of tip labels beneath it. The collection of these sets is a relabeling-
    invariant signature of the topology.
    """
    pos = 0
    clades: set[frozenset[str]] = set()

    def parse() -> set[str]:
        nonlocal pos
        tips: set[str] = set()
        if newick[pos] == "(":
            pos += 1  # consume '('
            while True:
                child = parse()
                tips |= child
                if newick[pos] == ",":
                    pos += 1
                    continue
                if newick[pos] == ")":
                    pos += 1
                    break
            # Skip a branch length on this internal node if present.
            _skip_length()
            clades.add(frozenset(tips))
        else:
            # A tip: read its label up to ':', ',' or ')'.
            m = re.match(r"[^():,;]+", newick[pos:])
            label = m.group(0)
            pos += len(label)
            _skip_length()
            tips.add(label.split(":")[0])
        return tips

    def _skip_length():
        nonlocal pos
        if pos < len(newick) and newick[pos] == ":":
            m = re.match(r":[-0-9.eE+]+", newick[pos:])
            if m:
                pos += len(m.group(0))

    parse()
    return frozenset(clades)


def topology_frequencies(newicks: list[str]) -> Counter:
    """Count how often each distinct topology appears across the draws."""
    sigs = [_clade_signature(nw) for nw in newicks]
    return Counter(sigs)


def topology_mixing_summary(newicks: list[str]) -> dict:
    """Summarize topology sampling for a single chain.

    Returns the number of distinct topologies, the fraction of consecutive draws
    where the topology changed (the switch rate), and the posterior probability of
    the most frequent topology.
    """
    if not newicks:
        return {"distinct": 0, "switch_rate": 0.0, "map_prob": 0.0}
    sigs = [_clade_signature(nw) for nw in newicks]
    switches = sum(1 for a, b in zip(sigs[:-1], sigs[1:]) if a != b)
    freqs = Counter(sigs)
    map_count = max(freqs.values())
    return {
        "distinct": len(freqs),
        "switch_rate": switches / max(len(sigs) - 1, 1),
        "map_prob": map_count / len(sigs),
    }


def clade_posteriors(newicks: list[str]) -> dict[frozenset[str], float]:
    """Posterior probability of each clade across the draws."""
    n = len(newicks)
    if n == 0:
        return {}
    counts: Counter = Counter()
    for nw in newicks:
        for clade in _clade_signature(nw):
            counts[clade] += 1
    return {clade: c / n for clade, c in counts.items()}


def clade_agreement(newicks_a: list[str], newicks_b: list[str]) -> float:
    """Agreement between two chains' clade posteriors.

    Returns one minus the mean absolute difference in clade probability over the
    union of clades seen in either chain. A value near 1 means the chains agree on
    the posterior tree distribution; a low value flags a mixing or convergence
    problem that scalar Rhat might miss.
    """
    pa = clade_posteriors(newicks_a)
    pb = clade_posteriors(newicks_b)
    clades = set(pa) | set(pb)
    if not clades:
        return 1.0
    diff = sum(abs(pa.get(c, 0.0) - pb.get(c, 0.0)) for c in clades) / len(clades)
    return 1.0 - diff
