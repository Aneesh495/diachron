// Lexical matrix view: languages (rows) by characters (columns), each cell
// shaded present / absent / missing. Missing is drawn distinctly from absent
// (a hatch, not a blank) so the viewer sees where there is no data versus where
// a cognate class is genuinely absent. Selecting languages in another view
// highlights their rows here, closing the linkage loop.

import type { Store } from "../store";

export interface MatrixData {
  languageLabels: string[];
  // columns[c][langIndex] in {1 present, 0 absent, -1 missing}
  columns: number[][];
}

export class MatrixView {
  private canvas: HTMLCanvasElement;
  constructor(
    private store: Store,
    container: HTMLElement,
    private data: MatrixData
  ) {
    this.canvas = document.createElement("canvas");
    this.canvas.className = "matrix-view";
    container.appendChild(this.canvas);
    this.store.subscribe(() => this.render());
    this.render();
  }

  private render(): void {
    const nLang = this.data.languageLabels.length;
    const nChar = Math.min(this.data.columns.length, 400); // cap for display
    const cellW = 3;
    const cellH = 8;
    this.canvas.width = nChar * cellW + 120;
    this.canvas.height = nLang * cellH + 20;
    const ctx = this.canvas.getContext("2d");
    if (!ctx) return;
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    const sel = this.store.selection;
    for (let c = 0; c < nChar; c += 1) {
      const col = this.data.columns[c];
      for (let l = 0; l < nLang; l += 1) {
        const state = col[l];
        const x = 120 + c * cellW;
        const y = 10 + l * cellH;
        if (state === 1) {
          ctx.fillStyle = sel.taxa.has(l) ? "#d6604d" : "#2c6e9e";
          ctx.fillRect(x, y, cellW, cellH - 1);
        } else if (state === 0) {
          ctx.fillStyle = "#e8eef3";
          ctx.fillRect(x, y, cellW, cellH - 1);
        } else {
          // missing: diagonal hatch so it is clearly not "absent"
          ctx.strokeStyle = "#cccccc";
          ctx.lineWidth = 0.5;
          ctx.beginPath();
          ctx.moveTo(x, y + cellH - 1);
          ctx.lineTo(x + cellW, y);
          ctx.stroke();
        }
      }
    }

    // Row labels, highlighting selected languages.
    ctx.font = "8px sans-serif";
    ctx.textBaseline = "middle";
    for (let l = 0; l < nLang; l += 1) {
      ctx.fillStyle = sel.taxa.has(l) ? "#d6604d" : "#333";
      ctx.fillText(this.data.languageLabels[l].slice(0, 18), 2, 10 + l * cellH + cellH / 2);
    }

    // Legend.
    ctx.font = "9px sans-serif";
    ctx.fillStyle = "#333";
    ctx.fillText("present / absent / missing (hatched)", 120, this.canvas.height - 4);
  }
}
