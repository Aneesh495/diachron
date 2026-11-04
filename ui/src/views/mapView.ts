// Map marks come from source records or exported posterior draws.
import type { Store } from "../store";
import type { AttestationPoint, GeoPoint, RootLocation } from "../types";
import { equirectangular, orthographic, longitude180, type Projection } from "../geography";

export class MapView {
  private el: SVGSVGElement;
  private globe = false;
  private centerLat = 35;
  private centerLon = 25;
  constructor(
    private store: Store,
    container: HTMLElement,
    private modernPoints: GeoPoint[],
    private attestations: AttestationPoint[]
  ) {
    const controls = document.createElement("div");
    controls.className = "map-controls";
    const projection = document.createElement("select");
    projection.setAttribute("aria-label", "Geographic projection");
    for (const [id, label] of [["map", "World map"], ["globe", "Globe"]]) {
      const option = document.createElement("option");
      option.value = id;
      option.textContent = label;
      projection.appendChild(option);
    }
    const rotations = document.createElement("span");
    rotations.hidden = true;
    for (const [label, caption, dlat, dlon] of [
      ["Rotate west", "←", 0, -30], ["Rotate east", "→", 0, 30],
      ["Rotate north", "↑", 20, 0], ["Rotate south", "↓", -20, 0],
    ] as const) {
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = caption;
      button.setAttribute("aria-label", label);
      button.addEventListener("click", () => {
        this.centerLat = Math.max(-80, Math.min(80, this.centerLat + dlat));
        this.centerLon = longitude180(this.centerLon + dlon);
        this.render();
      });
      rotations.appendChild(button);
    }
    projection.addEventListener("change", () => {
      this.globe = projection.value === "globe";
      rotations.hidden = !this.globe;
      this.render();
    });
    controls.append(projection, rotations);
    container.appendChild(controls);
    this.el = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    this.el.setAttribute("width", "100%");
    this.el.setAttribute("class", "map-view");
    this.el.setAttribute("role", "img");
    this.el.setAttribute("aria-label", "Geographic observations and sampled root locations");
    container.appendChild(this.el);
    this.store.subscribe(() => this.render());
    this.render();
  }

  private render(): void {
    while (this.el.firstChild) this.el.removeChild(this.el.firstChild);
    const w = this.el.parentElement?.clientWidth || 600;
    const h = Math.max(180, (this.el.parentElement?.clientHeight || 400) - 34);
    this.el.setAttribute("height", String(h));
    const globe = this.globe ? orthographic(w, h, this.centerLat, this.centerLon) : null;
    const project = globe?.project ?? equirectangular(w, h);
    const sel = this.store.selection;

    if (globe) {
      const sphere = this.svg("circle");
      sphere.setAttribute("cx", String(w / 2));
      sphere.setAttribute("cy", String(h / 2));
      sphere.setAttribute("r", String(globe.radius));
      sphere.setAttribute("fill", "#f5f8fa");
      sphere.setAttribute("stroke", "#b8c7d3");
      sphere.setAttribute("data-globe", "visible");
      this.el.appendChild(sphere);
    }
    for (let lon = -180; lon <= 150; lon += 30) {
      this.gridline(Array.from({ length: 61 }, (_, i) => project(i * 3 - 90, lon)));
    }
    for (let lat = -60; lat <= 60; lat += 30) {
      this.gridline(Array.from({ length: 121 }, (_, i) => project(lat, i * 3 - 180)));
    }

    for (const p of this.modernPoints) {
      const location = project(p.lat, p.lon);
      if (!location) continue;
      const [x, y] = location;
      const rect = this.svg("rect");
      rect.setAttribute("x", String(x - 3));
      rect.setAttribute("y", String(y - 3));
      rect.setAttribute("width", "6");
      rect.setAttribute("height", "6");
      rect.setAttribute("fill", "none");
      const taxon = this.store.bundle.tip_labels.indexOf(p.label);
      rect.setAttribute("stroke", sel.taxa.has(taxon) ? "#182739" : "#888");
      rect.setAttribute("stroke-width", "1");
      const title = this.svg("title");
      title.textContent = `modern reference: ${p.label} (present-day coordinate)`;
      rect.appendChild(title);
      this.el.appendChild(rect);
    }

    // An interval marks the years in which the observation could have occurred.
    for (const a of this.attestations) {
      if (a.lat === null || a.lon === null || sel.timeYbp < a.ageLow || sel.timeYbp > a.ageHigh) continue;
      const taxon = this.store.bundle.tip_labels.indexOf(a.label);
      if (sel.taxa.size > 0 && !sel.taxa.has(taxon)) continue;
      const location = project(a.lat, a.lon);
      if (!location) continue;
      const [x, y] = location;
      if (a.uncertaintyKm) {
        const ring = this.svg("circle");
        ring.setAttribute("cx", String(x));
        ring.setAttribute("cy", String(y));
        ring.setAttribute("r", String(Math.min(20, a.uncertaintyKm / 20)));
        ring.setAttribute("fill", "#e6a23c");
        ring.setAttribute("fill-opacity", "0.12");
        ring.setAttribute("stroke", "#e6a23c");
        ring.setAttribute("stroke-opacity", "0.4");
        this.el.appendChild(ring);
      }
      const dia = this.svg("path");
      dia.setAttribute("d", `M ${x} ${y - 4} L ${x + 4} ${y} L ${x} ${y + 4} L ${x - 4} ${y} Z`);
      dia.setAttribute("fill", "#e6a23c");
      const t = this.svg("title");
      t.textContent = `attestation: ${a.label}, ${a.ageLow}-${a.ageHigh} ybp (observed)`;
      dia.appendChild(t);
      this.el.appendChild(dia);
    }

    // Inferred latent cloud for the selected clade's ancestor at this time.
    this.renderLatentCloud(project);

    this.renderLegend(w);
  }

  private selectedLocations(): RootLocation[] {
    const { selection, bundle } = this.store;
    if (!selection.taxa.size || selection.taxa.size === bundle.tip_labels.length) {
      return bundle.rootLocations;
    }
    if (selection.clade) {
      const key = [...selection.clade].sort((a, b) => a - b).join(",");
      return bundle.cladeLocations[key] ?? [];
    }
    return [];
  }

  private renderLatentCloud(project: Projection): void {
    const selection = this.store.selection;
    for (const { lat, lon, age } of this.selectedLocations()) {
      if (Math.abs(age - selection.timeYbp) > 250 || Math.abs(lat) > 90 || Math.abs(lon) > 180) continue;
      const location = project(lat, lon);
      if (!location) continue;
      const [x, y] = location;
      const dot = this.svg("circle");
      dot.setAttribute("cx", String(x));
      dot.setAttribute("cy", String(y));
      dot.setAttribute("r", "3");
      dot.setAttribute("fill", "#4a6fa5");
      dot.setAttribute("fill-opacity", "0.35");
      this.el.appendChild(dot);
    }
  }

  private gridline(points: ReturnType<Projection>[]): void {
    let line = "";
    let drawing = false;
    for (const point of points) {
      if (point) line += `${drawing ? "L" : "M"} ${point[0]} ${point[1]} `;
      drawing = point !== null;
    }
    const grid = this.svg("path");
    grid.setAttribute("d", line);
    grid.setAttribute("fill", "none");
    grid.setAttribute("stroke", "#e8edf2");
    this.el.appendChild(grid);
  }

  private renderLegend(w: number): void {
    const items = [
      ["#888", "modern reference (present-day)"],
      ["#e6a23c", "attestation (observed)"],
      ["#4a6fa5", this.selectedLocations().length
        ? (this.store.selection.clade && this.store.selection.clade.size !== this.store.bundle.tip_labels.length
          ? "sampled clade (±250 years)" : "sampled root (±250 years)")
        : "location draws unavailable"],
    ];
    items.forEach(([color, label], i) => {
      const y = 16 + i * 16;
      const sw = this.svg("rect");
      sw.setAttribute("x", String(Math.max(4, w - 175)));
      sw.setAttribute("y", String(y - 8));
      sw.setAttribute("width", "10");
      sw.setAttribute("height", "10");
      sw.setAttribute("fill", color);
      this.el.appendChild(sw);
      const txt = this.svg("text");
      txt.setAttribute("x", String(Math.max(19, w - 160)));
      txt.setAttribute("y", String(y + 1));
      txt.setAttribute("font-size", "9");
      txt.setAttribute("fill", "#333");
      txt.textContent = label;
      this.el.appendChild(txt);
    });
  }

  private svg(tag: string): SVGElement {
    return document.createElementNS("http://www.w3.org/2000/svg", tag);
  }
}
