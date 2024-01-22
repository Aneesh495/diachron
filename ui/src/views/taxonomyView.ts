import type { Store } from "../store";
import type { TaxonomyGroup } from "../types";

export function taxonomyView(groups: TaxonomyGroup[], currentStore: () => Store): HTMLElement {
  const section = document.createElement("details");
  section.className = "taxonomy-comparison";
  const title = document.createElement("summary");
  title.textContent = "Compare catalog groups";
  const selector = document.createElement("select");
  selector.setAttribute("aria-label", "Catalog group");
  for (const group of groups) {
    const option = document.createElement("option");
    option.value = group.id;
    option.textContent = `${group.taxonomy} · ${group.label} (${group.taxa.length})`;
    selector.appendChild(option);
  }
  const details = document.createElement("div");
  section.append(title, selector, details);
  const render = () => {
    const group = groups.find((item) => item.id === selector.value);
    if (!group) return;
    details.replaceChildren();
    const source = document.createElement("a");
    source.textContent = `Source group in ${group.taxonomy}`;
    if (group.url.startsWith("https://")) source.href = group.url;
    source.target = "_blank";
    source.rel = "noopener noreferrer";
    const label = document.createElement("p");
    label.textContent = "Frequency of this exact group among the observed study tips.";
    const rows = document.createElement("ul");
    for (const [name, draws] of Object.entries(group.models)) {
      const row = document.createElement("li");
      const frequency = draws.clade_draws / draws.total_draws;
      row.textContent = `${name} · ${(frequency * 100).toFixed(1)}% · ${draws.clade_draws} of ${draws.total_draws} draws`;
      rows.appendChild(row);
    }
    const select = document.createElement("button");
    select.textContent = "Select observed members";
    select.addEventListener("click", () => currentStore().selectClade(new Set(group.taxa)));
    details.append(source, label, rows, select);
  };
  selector.addEventListener("change", render);
  render();
  return section;
}
