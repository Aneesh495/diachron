// Map marks come from source records or exported posterior draws.
import type { Store } from "../store";
import type { AttestationPoint, GeoPoint } from "../types";

export class MapView {
  private el: SVGSVGElement;
  constructor(
    private store: Store,
    container: HTMLElement,
    private modernPoints: GeoPoint[],
    private attestations: AttestationPoint[]
  ) {
    this.el = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    this.el.setAttribute("width", "100%");
    this.el.setAttribute("height", "100%");
    this.el.setAttribute("class", "map-view");
    container.appendChild(this.el);
    this.store.subscribe(() => this.render());
    this.render();
  }

  private project(lat: number, lon: number, w: number, h: number): [number, number] {
    return [((lon + 180) / 360) * w, ((90 - lat) / 180) * h];
  }

  private render(): void {
    while (this.el.firstChild) this.el.removeChild(this.el.firstChild);
    const w = this.el.clientWidth || 600;
    const h = this.el.clientHeight || 400;
    const sel = this.store.selection;

    for (const lon of [-120, -60, 0, 60, 120]) {
      const [x] = this.project(0, lon, w, h);
      this.gridline(x, 0, x, h);
    }
    for (const lat of [-60, -30, 0, 30, 60]) {
      const [, y] = this.project(lat, 0, w, h);
      this.gridline(0, y, w, y);
    }

    for (const p of this.modernPoints) {
      const [x, y] = this.project(p.lat, p.lon, w, h);
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
      const [x, y] = this.project(a.lat, a.lon, w, h);
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
    this.renderLatentCloud(w, h);

    this.renderLegend(w);
  }

  private renderLatentCloud(w: number, h: number): void {
    const selection = this.store.selection;
    if (selection.taxa.size > 0 && selection.taxa.size !== this.store.bundle.tip_labels.length) return;
    for (const { lat, lon, age } of this.store.bundle.rootLocations) {
      if (Math.abs(age - selection.timeYbp) > 250 || Math.abs(lat) > 90 || Math.abs(lon) > 180) continue;
      const [x, y] = this.project(lat, lon, w, h);
      const dot = this.svg("circle");
      dot.setAttribute("cx", String(x));
      dot.setAttribute("cy", String(y));
      dot.setAttribute("r", "3");
      dot.setAttribute("fill", "#4a6fa5");
      dot.setAttribute("fill-opacity", "0.35");
      this.el.appendChild(dot);
    }
  }

  private gridline(x1: number, y1: number, x2: number, y2: number): void {
    const line = this.svg("line");
    line.setAttribute("x1", String(x1));
    line.setAttribute("x2", String(x2));
    line.setAttribute("y1", String(y1));
    line.setAttribute("y2", String(y2));
    line.setAttribute("stroke", "#e8edf2");
    this.el.appendChild(line);
  }

  private renderLegend(w: number): void {
    const items = [
      ["#888", "modern reference (present-day)"],
      ["#e6a23c", "attestation (observed)"],
      ["#4a6fa5", this.store.bundle.rootLocations.length
        ? "sampled root (±250 years)"
        : "root location unavailable"],
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
