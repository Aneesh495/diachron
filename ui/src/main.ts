// Diachron workbench entry point.
//
// Loads a posterior bundle plus the dataset's geographic and matrix context,
// builds the shared store, and wires the four linked views (tree, map, timeline,
// matrix). A model switcher lets the viewer compare tree-only against
// borrowing-aware results side by side. Everything reads saved artifacts; the
// workbench never starts a chain.

import type { PosteriorBundle, GeoPoint } from "./types";
import { Store } from "./store";
import { TreeView } from "./views/treeView";
import { MapView } from "./views/mapView";
import { TimelineView } from "./views/timelineView";
import { MatrixView, type MatrixData } from "./views/matrixView";

interface WorkbenchData {
  bundles: Record<string, PosteriorBundle>;
  modernPoints: GeoPoint[];
  attestations: GeoPoint[];
  matrix: MatrixData;
}

export async function loadWorkbench(dataUrl: string): Promise<WorkbenchData> {
  const res = await fetch(dataUrl);
  if (!res.ok) throw new Error(`failed to load workbench data: ${res.status}`);
  return (await res.json()) as WorkbenchData;
}

export function mount(root: HTMLElement, data: WorkbenchData): void {
  const modelNames = Object.keys(data.bundles);
  let current = modelNames[0];
  let store = new Store(data.bundles[current]);

  root.innerHTML = "";
  root.className = "workbench";

  // Header with the honesty banner and model switcher.
  const header = document.createElement("div");
  header.className = "wb-header";
  const title = document.createElement("span");
  title.textContent = "Diachron workbench";
  title.className = "wb-title";
  header.appendChild(title);

  const syntheticFlag = document.createElement("span");
  syntheticFlag.className = "wb-synthetic";
  syntheticFlag.textContent = store.bundle.is_synthetic
    ? "synthetic known-truth corpus (not observed history)"
    : "observed-data study";
  header.appendChild(syntheticFlag);

  const select = document.createElement("select");
  for (const name of modelNames) {
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = `model: ${name}`;
    select.appendChild(opt);
  }
  header.appendChild(select);
  root.appendChild(header);

  const grid = document.createElement("div");
  grid.className = "wb-grid";
  root.appendChild(grid);

  const panels = ["tree", "map", "timeline", "matrix"].map((id) => {
    const p = document.createElement("div");
    p.className = `wb-panel wb-${id}`;
    const h = document.createElement("div");
    h.className = "wb-panel-title";
    h.textContent = id;
    p.appendChild(h);
    const body = document.createElement("div");
    body.className = "wb-panel-body";
    p.appendChild(body);
    grid.appendChild(p);
    return body;
  });

  const build = () => {
    for (const p of panels) p.innerHTML = "";
    new TreeView(store, panels[0]);
    new MapView(store, panels[1], data.modernPoints, data.attestations);
    new TimelineView(store, panels[2]);
    new MatrixView(store, panels[3], data.matrix);
  };
  build();

  select.addEventListener("change", () => {
    current = select.value;
    store = new Store(data.bundles[current]);
    syntheticFlag.textContent = store.bundle.is_synthetic
      ? "synthetic known-truth corpus (not observed history)"
      : "observed-data study";
    build();
  });
}

// Auto-mount when a root element and data URL are present.
const rootEl = document.getElementById("app");
if (rootEl) {
  const url = rootEl.dataset.bundle ?? "workbench.json";
  loadWorkbench(url)
    .then((data) => mount(rootEl, data))
    .catch((err) => {
      rootEl.textContent = `Could not load workbench data from ${url}: ${err}`;
    });
}
