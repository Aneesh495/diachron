import type { PosteriorBundle } from "../types";

export function modelComparison(bundles: Record<string, PosteriorBundle>): HTMLElement {
  const section = document.createElement("details");
  section.className = "model-comparison";
  const summary = document.createElement("summary");
  summary.textContent = "Compare saved models";
  const table = document.createElement("table");
  const title = table.createCaption();
  title.textContent = "Root ages are central intervals of displayed draws, in years before present.";
  const header = table.createTHead().insertRow();
  for (const text of ["Model", "Root age", "Root split R̂", "Bulk ESS", "Status"]) {
    const cell = document.createElement("th");
    cell.scope = "col";
    cell.textContent = text;
    header.appendChild(cell);
  }
  const body = table.createTBody();
  for (const [name, bundle] of Object.entries(bundles)) {
    const ages = bundle.chains.flatMap((chain) => chain.map((draw) => draw.root_age)).sort((a, b) => a - b);
    const quantile = (p: number) => ages[Math.floor(p * (ages.length - 1))].toFixed(0);
    const { split_rhat, bulk_ess } = bundle.rootAgeDiagnostic;
    const row = body.insertRow();
    const label = document.createElement("th");
    label.scope = "row";
    label.textContent = name;
    row.appendChild(label);
    for (const text of [
      `${quantile(0.5)} [${quantile(0.05)}, ${quantile(0.95)}]`,
      split_rhat?.toFixed(2) ?? "Unavailable",
      bulk_ess.toFixed(0),
      bundle.converged ? "Recorded scalar checks passed" : "Diagnostic draws",
    ]) row.insertCell().textContent = text;
  }
  section.append(summary, table);
  return section;
}
