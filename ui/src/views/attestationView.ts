import type { Store } from "../store";
import type { AttestationPoint } from "../types";

export class AttestationView {
  constructor(store: Store, container: HTMLElement, attestations: AttestationPoint[]) {
    const section = document.createElement("details");
    section.className = "attestation-inspector";
    const heading = document.createElement("summary");
    heading.textContent = "Attestation inspector";
    const select = document.createElement("select");
    select.setAttribute("aria-label", "Select an attestation");
    const record = document.createElement("div");
    section.append(heading, select, record);
    container.appendChild(section);

    const renderRecord = () => {
      record.replaceChildren();
      const item = attestations.find((a) => a.id === select.value);
      if (!item) {
        record.textContent = "No attestations for the current selection.";
        return;
      }
      const details = document.createElement("p");
      details.textContent = [
        `${item.region || "Unspecified locality"}, ${item.ageLow}–${item.ageHigh} years BP.`,
        item.locationNote || "", item.languageNote || "",
        item.conflicts.length ? `Alternatives: ${item.conflicts.join(", ")}.` : "",
      ].filter(Boolean).join(" ");
      record.appendChild(details);
      const source = document.createElement(item.sourceRef.startsWith("https://") ? "a" : "span");
      source.textContent = item.sourceRef;
      if (source instanceof HTMLAnchorElement) {
        source.href = item.sourceRef;
        source.rel = "noopener noreferrer";
        source.target = "_blank";
      }
      record.appendChild(source);
      for (const citation of [item.page, ...item.bibliography]) {
        if (!citation) continue;
        const row = document.createElement("p");
        row.textContent = citation;
        record.appendChild(row);
      }
      const jump = document.createElement("button");
      jump.type = "button";
      jump.textContent = "Show interval on map";
      jump.addEventListener("click", () => {
        store.setTime((item.ageLow + item.ageHigh) / 2);
        select.focus();
      });
      record.appendChild(jump);
    };

    const renderOptions = () => {
      const previous = select.value;
      const taxa = store.selection.taxa;
      const filtered = attestations.filter((a) =>
        taxa.size === 0 || taxa.has(store.bundle.tip_labels.indexOf(a.label))
      );
      select.replaceChildren(...filtered.map((item) => {
        const option = document.createElement("option");
        option.value = item.id;
        option.textContent = `${item.id} · ${item.label}`;
        return option;
      }));
      if (filtered.some((a) => a.id === previous)) select.value = previous;
      select.disabled = filtered.length === 0;
      renderRecord();
    };
    select.addEventListener("change", renderRecord);
    store.subscribe(renderOptions);
    renderOptions();
  }
}
