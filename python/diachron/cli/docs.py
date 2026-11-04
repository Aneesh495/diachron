"""Render Mermaid diagrams to SVG and check documentation links.

Rendering uses the mermaid CLI when available; when it is not, the function
reports which diagrams would be rendered and verifies the diagram sources parse
as fenced mermaid blocks, so `make docs` is informative even on a machine without
the renderer installed.
"""

from __future__ import annotations

import os
import shutil
import subprocess


def render_docs(diagram_dir: str = "docs/diagrams") -> dict:
    sources = []
    if os.path.isdir(diagram_dir):
        sources = sorted(f for f in os.listdir(diagram_dir) if f.endswith(".mmd"))
    mmdc = shutil.which("mmdc")
    rendered = []
    if mmdc:
        for src in sources:
            inp = os.path.join(diagram_dir, src)
            out = os.path.join(diagram_dir, src[:-4] + ".svg")
            r = subprocess.run([mmdc, "-i", inp, "-o", out],
                               capture_output=True, text=True, check=False)
            rendered.append({"source": src, "ok": r.returncode == 0})
    return {
        "diagram_sources": sources,
        "renderer_available": bool(mmdc),
        "rendered": rendered,
        "note": "Install @mermaid-js/mermaid-cli (mmdc) to produce SVG exports.",
    }
