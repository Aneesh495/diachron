// Tree view: renders the maximum-credibility tree as an SVG cladogram with
// branch lengths in years, clade support annotations, and clickable clades.
//
// Internal nodes are drawn with their posterior support so the viewer never
// mistakes a weakly supported split for a certain one; low-support nodes are
// drawn faintly. Clicking a node selects its clade and propagates to the other
// views through the store.

import type { Store } from "../store";
import type { ParsedTree } from "../types";
import { descendantTaxa } from "../newick";

export class TreeView {
  private el: SVGSVGElement;
  private tree: ParsedTree;
  private fitted: boolean = true;
  constructor(private store: Store, container: HTMLElement) {
    this.tree = store.summary.maxCredibilityTree();
    const controls = document.createElement("div");
    controls.className = "tree-controls";
    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.textContent = "Show tip labels";
    toggle.setAttribute("aria-pressed", "false");
    toggle.addEventListener("click", () => {
      this.fitted = !this.fitted;
      toggle.textContent = this.fitted ? "Show tip labels" : "Fit full tree";
      toggle.setAttribute("aria-pressed", String(!this.fitted));
      container.scrollTop = 0;
      this.render();
    });
    const select = document.createElement("button");
    select.type = "button";
    select.textContent = "Select family";
    select.addEventListener("click", () => store.selectClade(descendantTaxa(this.tree, this.tree.root)));
    controls.append(toggle, select);
    container.appendChild(controls);

    this.el = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    this.el.setAttribute("width", "100%");
    this.el.setAttribute("class", "tree-view");
    container.appendChild(this.el);
    this.store.subscribe(() => this.render());
    this.render();
  }

  private render(): void {
    const tree = this.tree;
    while (this.el.firstChild) this.el.removeChild(this.el.firstChild);

    const width = this.el.parentElement?.clientWidth || 600;
    const viewport = Math.max(220, (this.el.parentElement?.clientHeight ?? 300) - 34);
    const margin = 24;
    const rootAge = tree.nodes[tree.root].height || 1;

    // Assign a y position to each tip in traversal order, x by node age.
    const tipOrder: number[] = [];
    const collectTips = (node: number) => {
      const kids = tree.nodes[node].children;
      if (kids.length === 0) {
        tipOrder.push(node);
      } else {
        for (const c of kids) collectTips(c);
      }
    };
    collectTips(tree.root);

    const height = this.fitted ? viewport : Math.max(viewport, tipOrder.length * 16 + 2 * margin);
    this.el.setAttribute("height", String(height));
    const showLabels = !this.fitted || tipOrder.length * 12 < height;
    const right = showLabels ? Math.min(width / 2, 175) : margin;

    const yOf = new Map<number, number>();
    tipOrder.forEach((n, i) => {
      yOf.set(n, margin + (i * (height - 2 * margin)) / Math.max(tipOrder.length - 1, 1));
    });
    // Internal node y is the mean of its children.
    const order = this.postorder(tree);
    for (const node of order) {
      const kids = tree.nodes[node].children;
      if (kids.length) {
        const mean = kids.reduce((s, c) => s + (yOf.get(c) ?? 0), 0) / kids.length;
        yOf.set(node, mean);
      }
    }

    const xOf = (age: number) =>
      width - right - (age / rootAge) * (width - right - margin);

    const sel = this.store.selection;

    // Draw edges.
    for (const node of tree.nodes) {
      if (node.parent === null) continue;
      const p = tree.nodes[node.parent];
      const x1 = xOf(p.height);
      const x2 = xOf(node.height);
      const y1 = yOf.get(node.id)!;
      const yParent = yOf.get(p.id)!;
      const inClade =
        node.taxon !== null && sel.taxa.has(node.taxon);
      const line = this.svg("path");
      line.setAttribute(
        "d",
        `M ${x1} ${yParent} L ${x1} ${y1} L ${x2} ${y1}`
      );
      line.setAttribute("fill", "none");
      line.setAttribute("stroke", inClade ? "#d6604d" : "#4a6fa5");
      line.setAttribute("stroke-width", inClade ? "2.5" : "1.4");
      this.el.appendChild(line);
    }

    // Draw internal-node support dots and tip labels.
    for (const node of tree.nodes) {
      const x = xOf(node.height);
      const y = yOf.get(node.id)!;
      if (node.children.length) {
        const taxa = descendantTaxa(tree, node.id);
        const support = this.store.summary.supportForTaxa(taxa);
        const dot = this.svg("circle");
        dot.setAttribute("cx", String(x));
        dot.setAttribute("cy", String(y));
        dot.setAttribute("r", "4");
        // Support drives opacity: weak clades are faint, never shown as certain.
        dot.setAttribute("fill", "#333");
        dot.setAttribute("fill-opacity", String(0.15 + 0.85 * support));
        dot.style.cursor = "pointer";
        dot.addEventListener("click", () => this.store.selectClade(taxa));
        const title = this.svg("title");
        title.textContent = `clade frequency ${(support * 100).toFixed(0)}%, age ${node.height.toFixed(0)} ybp`;
        dot.appendChild(title);
        this.el.appendChild(dot);
      } else if (showLabels) {
        const label = this.svg("text");
        label.setAttribute("x", String(x + 6));
        label.setAttribute("y", String(y + 3));
        label.setAttribute("font-size", "11");
        label.setAttribute(
          "fill",
          sel.taxa.has(node.taxon!) ? "#d6604d" : "#222"
        );
        label.textContent = this.store.bundle.tip_names[node.label ?? ""] ?? node.label ?? `t${node.taxon}`;
        label.style.cursor = "pointer";
        label.addEventListener("click", () => this.store.toggleTaxon(node.taxon!));
        this.el.appendChild(label);
      }
    }

    // Axis caption with the synthetic-data honesty note.
    const caption = this.svg("text");
    caption.setAttribute("x", String(margin));
    caption.setAttribute("y", String(height - 8));
    caption.setAttribute("font-size", "10");
    caption.setAttribute("fill", "#666");
    caption.textContent = `max-credibility tree, branch lengths in years${
      this.store.bundle.is_synthetic ? " (synthetic known-truth corpus)" : ""
    }`;
    this.el.appendChild(caption);
  }

  private postorder(tree: ParsedTree): number[] {
    const out: number[] = [];
    const stack: Array<[number, boolean]> = [[tree.root, false]];
    while (stack.length) {
      const [node, done] = stack.pop()!;
      if (done) out.push(node);
      else {
        stack.push([node, true]);
        for (const c of tree.nodes[node].children) stack.push([c, false]);
      }
    }
    return out;
  }

  private svg(tag: string): SVGElement {
    return document.createElementNS("http://www.w3.org/2000/svg", tag);
  }
}
