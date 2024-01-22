// Posterior summaries computed in the browser from a loaded bundle.
//
// These read the already-sampled draws. Clade support is the posterior
// probability of each bipartition across draws; the consensus tree keeps clades
// above a support threshold. Scalar summaries give credible intervals for the
// timeline and legends. Nothing here samples a new chain.

import type { PosteriorBundle, ParsedTree } from "./types";
import { parseNewick, descendantTaxa } from "./newick";

export interface CladeSupport {
  taxa: Set<number>;
  support: number; // posterior probability in [0, 1]
}

export interface ScalarSummary {
  mean: number;
  sd: number;
  q05: number;
  q50: number;
  q95: number;
}

export class PosteriorSummary {
  readonly bundle: PosteriorBundle;
  readonly trees: ParsedTree[];
  readonly cladeSupport: CladeSupport[];

  constructor(bundle: PosteriorBundle) {
    this.bundle = bundle;
    this.trees = [];
    for (const chain of bundle.chains) {
      for (const draw of chain) {
        this.trees.push(parseNewick(draw.newick, bundle.tip_labels, draw.root_age));
      }
    }
    this.cladeSupport = this.computeCladeSupport();
  }

  private computeCladeSupport(): CladeSupport[] {
    const counts = new Map<string, { taxa: Set<number>; n: number }>();
    for (const tree of this.trees) {
      for (const node of tree.nodes) {
        if (node.children.length === 0) continue; // tips are trivial
        const taxa = descendantTaxa(tree, node.id);
        if (taxa.size < 2) continue;
        const key = [...taxa].sort((a, b) => a - b).join(",");
        const entry = counts.get(key);
        if (entry) entry.n += 1;
        else counts.set(key, { taxa, n: 1 });
      }
    }
    const total = this.trees.length || 1;
    return [...counts.values()]
      .map((e) => ({ taxa: e.taxa, support: e.n / total }))
      .sort((a, b) => b.support - a.support);
  }

  // Credible interval and moments for a scalar draw field.
  scalarSummary(field: "root_age" | "gain" | "loss" | "clock_rate"): ScalarSummary {
    const values: number[] = [];
    for (const chain of this.bundle.chains) {
      for (const d of chain) values.push(d[field]);
    }
    values.sort((a, b) => a - b);
    const n = values.length;
    const mean = values.reduce((a, b) => a + b, 0) / n;
    const variance = values.reduce((a, b) => a + (b - mean) ** 2, 0) / Math.max(n - 1, 1);
    const q = (p: number) => values[Math.min(n - 1, Math.floor(p * n))];
    return { mean, sd: Math.sqrt(variance), q05: q(0.05), q50: q(0.5), q95: q(0.95) };
  }

  // The highest-support tree, used as the display consensus. Returns the sampled
  // tree whose clades have the greatest summed support, a maximum-credibility
  // choice that is guaranteed to be a real sampled tree (not a synthetic merge).
  maxCredibilityTree(): ParsedTree {
    const supportByKey = new Map<string, number>();
    for (const cs of this.cladeSupport) {
      const key = [...cs.taxa].sort((a, b) => a - b).join(",");
      supportByKey.set(key, cs.support);
    }
    let best = this.trees[0];
    let bestScore = Number.NEGATIVE_INFINITY;
    for (const tree of this.trees) {
      let score = 0;
      for (const node of tree.nodes) {
        if (node.children.length === 0) continue;
        const taxa = descendantTaxa(tree, node.id);
        const key = [...taxa].sort((a, b) => a - b).join(",");
        score += Math.log(Math.max(supportByKey.get(key) ?? 1e-6, 1e-6));
      }
      if (score > bestScore) {
        bestScore = score;
        best = tree;
      }
    }
    return best;
  }

  // Support for the clade defined by exactly a set of taxa (for the legend).
  supportForTaxa(taxa: Set<number>): number {
    const key = [...taxa].sort((a, b) => a - b).join(",");
    for (const cs of this.cladeSupport) {
      const csKey = [...cs.taxa].sort((a, b) => a - b).join(",");
      if (csKey === key) return cs.support;
    }
    return 0;
  }
}
