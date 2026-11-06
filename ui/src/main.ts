// Diachron workbench entry point.
//
// Loads a posterior bundle plus the dataset's geographic and matrix context,
// builds the shared store, and wires the four linked views (tree, map, timeline,
// matrix). A model switcher lets the viewer compare tree-only against
// borrowing-aware results side by side. Everything reads saved artifacts; the
// workbench never starts a chain.

import type { AttestationPoint, PosteriorBundle, GeoPoint, TourStep, TaxonomyGroup } from "./types";
import { Store } from "./store";
import { modelComparison } from "./views/modelComparison";
import { AttestationView } from "./views/attestationView";
import { TreeView } from "./views/treeView";
import { MapView } from "./views/mapView";
import { TimelineView } from "./views/timelineView";
import { MatrixView, type MatrixData } from "./views/matrixView";
import { researchTour } from "./views/researchTour";
import { taxonomyView } from "./views/taxonomyView";

interface WorkbenchData {
  bundles: Record<string, PosteriorBundle>;
  modernPoints: GeoPoint[];
  attestations: AttestationPoint[];
  matrix: MatrixData;
  tour?: TourStep[];
  taxonomy?: TaxonomyGroup[];
}

export async function loadWorkbench(dataUrl: string): Promise<WorkbenchData> {
  const res = await fetch(dataUrl);
  if (!res.ok) throw new Error(`failed to load workbench data: ${res.status}`);
  return (await res.json()) as WorkbenchData;
}

export function mount(root: HTMLElement, data: WorkbenchData): void {
  const modelNames = Object.keys(data.bundles);
  if (!modelNames.length) throw new Error("No saved posterior models");
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
  const mixing = document.createElement("span");
  mixing.className = "wb-mixing";
  mixing.setAttribute("role", "status");
  const renderMixing = () => {
    const { converged, rootAgeDiagnostic } = store.bundle;
    if (converged) {
      mixing.textContent = "Recorded scalar checks passed";
      return;
    }
    const rhat = rootAgeDiagnostic.split_rhat?.toFixed(2) ?? "unavailable";
    mixing.textContent = `Diagnostic draws · root R̂ ${rhat} · bulk ESS ${rootAgeDiagnostic.bulk_ess.toFixed(0)}`;
  };
  renderMixing();
  header.appendChild(mixing);

  const select = document.createElement("select");
  select.setAttribute("aria-label", "Posterior model");
  for (const name of modelNames) {
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = `model: ${name}`;
    select.appendChild(opt);
  }
  header.appendChild(select);
  root.appendChild(header);

  const research = document.createElement("aside");
  research.className = "research-tools";
  root.appendChild(research);

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
    new AttestationView(store, panels[2], data.attestations);
    if (modelNames.length > 1) panels[2].appendChild(modelComparison(data.bundles));
    new MatrixView(store, panels[3], data.matrix);
  };
  build();

  const selectModel = (name: string) => {
    if (name === current || !modelNames.includes(name)) return;
    const previous = store.selection;
    const matrixFocus = store.matrixFocus;
    current = name;
    select.value = name;
    store = new Store(data.bundles[current]);
    store.selection = {
      taxa: new Set(previous.taxa), clade: previous.clade && new Set(previous.clade),
      timeYbp: previous.timeYbp,
    };
    store.matrixFocus = matrixFocus;
    syntheticFlag.textContent = store.bundle.is_synthetic
      ? "synthetic known-truth corpus (not observed history)"
      : "observed-data study";
    renderMixing();
    build();
  };
  select.addEventListener("change", () => selectModel(select.value));
  if (data.tour?.length) {
    const dated = modelNames.find((name) => data.bundles[name].model_name === "dated");
    research.appendChild(researchTour(data.tour, () => store, selectModel, dated));
  }
  if (data.taxonomy?.length) research.appendChild(taxonomyView(data.taxonomy, () => store));
  if (!research.childElementCount) research.remove();
}

// Auto-mount when a root element and data URL are present.
const rootEl = document.getElementById("app");
if (rootEl) {
  const url = rootEl.dataset.bundle ?? "workbench.json";
  const inline = document.getElementById("workbench-data");
  const data = inline
    ? Promise.resolve().then(() => JSON.parse(inline.textContent ?? "") as WorkbenchData)
    : loadWorkbench(url);
  data.then((snapshot) => mount(rootEl, snapshot))
    .catch((err) => {
      rootEl.textContent = `Could not load workbench data from ${url}: ${err}`;
    });
}
