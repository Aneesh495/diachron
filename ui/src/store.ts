// A tiny observable store holding the current selection and bundle.
//
// All four views subscribe to the store and re-render when the selection
// changes. Selecting a clade in the tree highlights the same languages on the
// globe and in the matrix, and moving the timeline filters what every view
// shows. The store is the single source of truth for cross-view linkage.

import type { PosteriorBundle, Selection } from "./types";
import { PosteriorSummary } from "./posterior";

type Listener = () => void;

export class Store {
  bundle: PosteriorBundle;
  summary: PosteriorSummary;
  selection: Selection;
  private listeners: Listener[] = [];

  constructor(bundle: PosteriorBundle) {
    this.bundle = bundle;
    this.summary = new PosteriorSummary(bundle);
    const rootAge = this.summary.scalarSummary("root_age");
    this.selection = {
      taxa: new Set(),
      clade: null,
      timeYbp: rootAge.q50,
    };
  }

  subscribe(fn: Listener): void {
    this.listeners.push(fn);
  }

  private emit(): void {
    for (const fn of this.listeners) fn();
  }

  selectClade(taxa: Set<number>): void {
    this.selection.clade = taxa;
    this.selection.taxa = new Set(taxa);
    this.emit();
  }

  toggleTaxon(taxon: number): void {
    if (this.selection.taxa.has(taxon)) this.selection.taxa.delete(taxon);
    else this.selection.taxa.add(taxon);
    this.selection.clade = null;
    this.emit();
  }

  clearSelection(): void {
    this.selection.taxa = new Set();
    this.selection.clade = null;
    this.emit();
  }

  setTime(ybp: number): void {
    this.selection.timeYbp = ybp;
    this.emit();
  }
}
