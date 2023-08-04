// Newick parser for the workbench.
//
// The bundle stores each sampled tree as a Newick string with branch lengths in
// years. Derive ages from the recorded root so historical terminal nodes
// retain their sampled ages.

import type { ParsedTree, TreeNode } from "./types";

export function parseNewick(newick: string, tipLabels: string[], rootAge: number): ParsedTree {
  const labelToTaxon = new Map<string, number>();
  tipLabels.forEach((lab, i) => labelToTaxon.set(lab, i));

  const nodes: TreeNode[] = [];
  let pos = 0;

  function addNode(taxon: number | null, branch: number, label?: string): number {
    const id = nodes.length;
    nodes.push({
      id,
      parent: null,
      children: [],
      branchLength: branch,
      height: 0,
      taxon,
      label,
    });
    return id;
  }

  function readLength(): number {
    if (pos < newick.length && newick[pos] === ":") {
      const m = /^:[-0-9.eE+]+/.exec(newick.slice(pos));
      if (m) {
        pos += m[0].length;
        return parseFloat(m[0].slice(1));
      }
    }
    return 0;
  }

  function parse(): number {
    if (newick[pos] === "(") {
      pos += 1;
      const kids: number[] = [];
      for (;;) {
        kids.push(parse());
        if (newick[pos] === ",") {
          pos += 1;
          continue;
        }
        if (newick[pos] === ")") {
          pos += 1;
          break;
        }
      }
      const branch = readLength();
      const node = addNode(null, branch);
      for (const k of kids) {
        nodes[k].parent = node;
        nodes[node].children.push(k);
      }
      return node;
    }
    const m = /^[^():,;]+/.exec(newick.slice(pos));
    const token = m ? m[0] : "";
    pos += token.length;
    const label = token.split(":")[0];
    let branch = 0;
    if (token.includes(":")) {
      branch = parseFloat(token.split(":")[1]);
    } else {
      branch = readLength();
    }
    const taxon = labelToTaxon.has(label) ? labelToTaxon.get(label)! : null;
    return addNode(taxon, branch, label);
  }

  const root = parse();

  nodes[root].height = rootAge;
  const stack = [root];
  while (stack.length) {
    const parent = stack.pop()!;
    for (const child of nodes[parent].children) {
      nodes[child].height = Math.max(0, nodes[parent].height - nodes[child].branchLength);
      stack.push(child);
    }
  }

  return { nodes, root };
}

// Taxa beneath a node, used to map a clade selection to tips.
export function descendantTaxa(tree: ParsedTree, node: number): Set<number> {
  const out = new Set<number>();
  const stack = [node];
  while (stack.length) {
    const n = stack.pop()!;
    const tx = tree.nodes[n].taxon;
    if (tx !== null) out.add(tx);
    for (const c of tree.nodes[n].children) stack.push(c);
  }
  return out;
}
