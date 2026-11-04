// Shared types for the Diachron workbench.
//
// The workbench reads posterior bundles produced by the Python inference layer.
// It never runs inference itself: moving the timeline or selecting a clade reads
// already-sampled draws, it does not trigger a new chain. These types mirror the
// bundle schema emitted by `diachron.infer.artifacts.PosteriorBundle`.

export interface Draw {
  iteration: number;
  log_posterior: number;
  root_age: number;
  gain: number;
  loss: number;
  clock_rate: number;
  newick: string;
}

export interface PosteriorBundle {
  model_name: string;
  dataset_name: string;
  dataset_version: string;
  is_synthetic: boolean;
  n_chains: number;
  draws_per_chain: number;
  tip_labels: string[];
  chains: Draw[][];
  source_commit: string;
  notes: string;
  rootLocations: RootLocation[];
}

// A parsed phylogenetic tree node. `provenance` records whether a node's
// position is observed, inferred, or illustrative, which drives the honesty
// legend: inferred internal nodes are never drawn as if they were observed.
export type Provenance = "observed" | "inferred" | "illustrative";

export interface TreeNode {
  id: number;
  parent: number | null;
  children: number[];
  branchLength: number;
  height: number; // years before present
  taxon: number | null; // tip taxon index, else null
  label?: string;
}

export interface ParsedTree {
  nodes: TreeNode[];
  root: number;
}

// A geographic point with an explicit kind so the map can keep modern reference
// points, historical attestations, and inferred latent locations visually
// distinct (gate G09).
export type GeoKind =
  | "modern-reference"
  | "attestation"
  | "inferred-latent"
  | "illustrative-spread";

export interface GeoPoint {
  lat: number;
  lon: number;
  kind: GeoKind;
  label: string;
  uncertaintyKm?: number;
}

export interface AttestationPoint {
  id: string;
  lat: number | null;
  lon: number | null;
  kind: "attestation";
  label: string;
  ageLow: number;
  ageHigh: number;
  uncertaintyKm: number | null;
  region: string | null;
  sourceRef: string;
  page: string | null;
  bibliography: string[];
  locationKind: string | null;
  locationNote: string | null;
  languageNote: string | null;
  conflicts: string[];
}

export interface RootLocation {
  lat: number;
  lon: number;
  age: number;
}

// The current selection, propagated across all linked views.
export interface Selection {
  taxa: Set<number>; // selected tip taxa
  clade: Set<number> | null; // a selected clade (set of taxa)
  timeYbp: number; // current timeline position, years before present
}
