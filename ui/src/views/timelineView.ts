// Timeline view: a slider over years before present with the posterior root-age
// credible interval drawn as a band. Moving the slider sets the store's time,
// which filters the map (attestations appear only once their date is reached)
// and annotates the tree. The slider reads sampled draws; it never triggers new
// inference, so dragging it is cheap and bounded.

import type { Store } from "../store";

export class TimelineView {
  private el: HTMLElement;
  constructor(private store: Store, container: HTMLElement) {
    this.el = document.createElement("div");
    this.el.className = "timeline-view";
    container.appendChild(this.el);
    this.store.subscribe(() => this.renderReadout());
    this.render();
  }

  private render(): void {
    const root = this.store.summary.scalarSummary("root_age");
    const maxAge = Math.ceil(root.q95 * 1.1);

    this.el.innerHTML = "";
    const band = document.createElement("div");
    band.className = "timeline-band";
    const lo = (1 - root.q95 / maxAge) * 100;
    const hi = (1 - root.q05 / maxAge) * 100;
    band.style.cssText =
      `position:relative;height:10px;border:1px solid #ccc;` +
      `background:linear-gradient(to right,` +
      `#eee 0%, #eee ${lo}%, #bcd ${lo}%, #bcd ${hi}%, #eee ${hi}%, #eee 100%);`;
    this.el.appendChild(band);

    const slider = document.createElement("input");
    slider.type = "range";
    slider.min = "0";
    slider.max = String(maxAge);
    slider.step = "50";
    slider.value = String(Math.round(this.store.selection.timeYbp));
    slider.style.width = "100%";
    slider.addEventListener("input", () => {
      this.store.setTime(Number(slider.value));
    });
    this.el.appendChild(slider);

    const readout = document.createElement("div");
    readout.className = "timeline-readout";
    readout.style.cssText = "font-size:11px;color:#333;margin-top:4px;";
    this.el.appendChild(readout);
    this.renderReadout();
  }

  private renderReadout(): void {
    const readout = this.el.querySelector(".timeline-readout");
    if (!readout) return;
    const root = this.store.summary.scalarSummary("root_age");
    readout.textContent =
      `time: ${Math.round(this.store.selection.timeYbp)} years before present  |  ` +
      `root age posterior: ${root.q50.toFixed(0)} ` +
      `[${root.q05.toFixed(0)}, ${root.q95.toFixed(0)}] (90% credible)`;
  }
}
