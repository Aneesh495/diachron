import type { Store } from "../store";
import type { TourStep } from "../types";

export function researchTour(
  steps: TourStep[],
  currentStore: () => Store,
  selectModel: (name: string) => void,
  dateModel: string | undefined,
): HTMLElement {
  const section = document.createElement("details");
  section.className = "research-tour";
  const heading = document.createElement("summary");
  heading.textContent = "Research tour";
  const navigation = document.createElement("div");
  navigation.className = "tour-navigation";
  const previous = document.createElement("button");
  previous.textContent = "Previous";
  const next = document.createElement("button");
  next.textContent = "Next";
  const position = document.createElement("span");
  position.setAttribute("aria-live", "polite");
  navigation.append(previous, position, next);
  const title = document.createElement("h2");
  const body = document.createElement("p");
  const evidence = document.createElement("div");
  section.append(heading, navigation, title, body, evidence);
  let index = 0;

  const link = (label: string, url: string): HTMLAnchorElement => {
    const anchor = document.createElement("a");
    anchor.textContent = label;
    if (url.startsWith("https://")) anchor.href = url;
    anchor.target = "_blank";
    anchor.rel = "noopener noreferrer";
    return anchor;
  };
  const render = () => {
    const step = steps[index];
    title.textContent = step.title;
    body.textContent = step.body;
    position.textContent = `${index + 1} of ${steps.length}`;
    previous.disabled = index === 0;
    next.disabled = index === steps.length - 1;
    evidence.replaceChildren();
    if (step.cells) {
      const cells = step.cells;
      const show = document.createElement("button");
      show.textContent = "Select example cells";
      show.addEventListener("click", () => currentStore().focusMatrix(cells));
      evidence.appendChild(show);
    }
    if (step.dates) {
      const taxon = step.taxon;
      for (const date of step.dates) {
        const row = document.createElement("p");
        row.textContent = `${date.label} · ${date.source_date}. `;
        row.appendChild(link(date.citation, date.url));
        const show = document.createElement("button");
        show.textContent = "Show approximate date";
        show.addEventListener("click", () => {
          if (dateModel) selectModel(dateModel);
          currentStore().selectClade(new Set(taxon === undefined ? [] : [taxon]), date.age);
        });
        row.append(" ", show);
        evidence.appendChild(row);
      }
    }
    if (step.model) {
      const model = step.model;
      const show = document.createElement("button");
      show.textContent = "Show root location draws";
      show.addEventListener("click", () => {
        selectModel(model);
        const store = currentStore();
        store.selectClade(new Set(store.bundle.tip_labels.map((_, index) => index)),
          store.summary.scalarSummary("root_age").q50);
      });
      evidence.appendChild(show);
    }
    for (const source of step.links) {
      const row = document.createElement("p");
      row.appendChild(link(source.label, source.url));
      evidence.appendChild(row);
    }
  };
  previous.addEventListener("click", () => { index -= 1; render(); });
  next.addEventListener("click", () => { index += 1; render(); });
  render();
  return section;
}
