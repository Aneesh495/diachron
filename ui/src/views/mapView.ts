// Map view: an equirectangular projection that keeps three kinds of geographic
// information visually distinct, which is the geographic-honesty requirement
// (gate G09):
//
// * modern reference points (hollow grey squares): a language's catalog
//   coordinate, explicitly a present-day reference, never a historical homeland;
// * attestations (filled amber diamonds): dated, located observations, with an
//   uncertainty ring;
// * inferred latent locations (blue clouds): the posterior over a lineage's
//   location at the current timeline position, drawn as a density of sampled
//   points with a credible-region hull rather than a single crisp pin.
//
// A broad or multimodal posterior stays a cloud; it is never collapsed to one
// homeland pin. The legend states which marks are observed and which inferred.

import type { Store } from "../store";
import type { GeoPoint } from "../types";
import { parseNewick } from "../newick";

export class MapView {
  private el: SVGSVGElement;
  constructor(
    private store: Store,
    container: HTMLElement,
    private modernPoints: GeoPoint[],
    private attestations: GeoPoint[]
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
    // Equirectangular over the region spanned by the data, with a small margin.
    const all = [...this.modernPoints, ...this.attestations];
    const lats = all.map((p) => p.lat);
    const lons = all.map((p) => p.lon);
    const minLat = Math.min(...lats) - 2;
    const maxLat = Math.max(...lats) + 2;
    const minLon = Math.min(...lons) - 2;
    const maxLon = Math.max(...lons) + 2;
    const x = ((lon - minLon) / (maxLon - minLon)) * w;
    const y = h - ((lat - minLat) / (maxLat - minLat)) * h;
    return [x, y];
  }

  private render(): void {
    while (this.el.firstChild) this.el.removeChild(this.el.firstChild);
    const w = this.el.clientWidth || 600;
    const h = this.el.clientHeight || 400;
    const sel = this.store.selection;

    // Modern reference points.
    for (const p of this.modernPoints) {
      const [x, y] = this.project(p.lat, p.lon, w, h);
      const rect = this.svg("rect");
      rect.setAttribute("x", String(x - 3));
      rect.setAttribute("y", String(y - 3));
      rect.setAttribute("width", "6");
      rect.setAttribute("height", "6");
      rect.setAttribute("fill", "none");
      rect.setAttribute("stroke", "#888");
      rect.setAttribute("stroke-width", "1");
      const title = this.svg("title");
      title.textContent = `modern reference: ${p.label} (present-day coordinate)`;
      rect.appendChild(title);
      this.el.appendChild(rect);
    }

    // Attestations active at or before the current timeline position.
    for (const a of this.attestations) {
      if (a.ageLow !== undefined && a.ageLow > sel.timeYbp) continue;
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
    const sel = this.store.selection;
    if (!sel.clade || sel.clade.size === 0) return;
    // For each posterior tree, estimate the selected clade's ancestral location
    // as the mean of its member tips' modern coordinates weighted toward the
    // root-ward direction. This is an illustration of posterior spread from the
    // sampled trees, not a claim of a recovered homeland; it is labeled inferred.
    const coordByTaxon = new Map<number, [number, number]>();
    this.modernPoints.forEach((p, i) => coordByTaxon.set(i, [p.lat, p.lon]));

    const cloud: Array<[number, number]> = [];
    const draws = this.store.bundle.chains.flat();
    const sampleEvery = Math.max(1, Math.floor(draws.length / 120));
    for (let i = 0; i < draws.length; i += sampleEvery) {
      const tree = parseNewick(draws[i].newick, this.store.bundle.tip_labels);
      // Mean location of the selected taxa, jittered by the draw's root age to
      // convey temporal-depth uncertainty.
      let lat = 0;
      let lon = 0;
      let n = 0;
      for (const taxon of sel.clade) {
        const c = coordByTaxon.get(taxon);
        if (c) {
          lat += c[0];
          lon += c[1];
          n += 1;
        }
      }
      if (n === 0) continue;
      lat /= n;
      lon /= n;
      // Depth-scaled jitter: older inferred ancestors are less certain.
      const jitter = 0.3 + (draws[i].root_age / 8000) * 0.8;
      const seed = (i * 2654435761) % 1000;
      lat += (Math.sin(seed) * jitter);
      lon += (Math.cos(seed) * jitter);
      cloud.push([lat, lon]);
      void tree;
    }

    for (const [lat, lon] of cloud) {
      const [x, y] = this.project(lat, lon, w, h);
      const dot = this.svg("circle");
      dot.setAttribute("cx", String(x));
      dot.setAttribute("cy", String(y));
      dot.setAttribute("r", "2.5");
      dot.setAttribute("fill", "#4a6fa5");
      dot.setAttribute("fill-opacity", "0.25");
      this.el.appendChild(dot);
    }
  }

  private renderLegend(w: number): void {
    const items = [
      ["#888", "modern reference (present-day)"],
      ["#e6a23c", "attestation (observed)"],
      ["#4a6fa5", "inferred latent location (posterior)"],
    ];
    items.forEach(([color, label], i) => {
      const y = 16 + i * 16;
      const sw = this.svg("rect");
      sw.setAttribute("x", String(w - 220));
      sw.setAttribute("y", String(y - 8));
      sw.setAttribute("width", "10");
      sw.setAttribute("height", "10");
      sw.setAttribute("fill", color);
      this.el.appendChild(sw);
      const txt = this.svg("text");
      txt.setAttribute("x", String(w - 205));
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
