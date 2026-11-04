"""Validate checked-in diagrams and optionally render fresh copies."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET


def _is_svg(path: Path) -> bool:
    try:
        return ET.parse(path).getroot().tag == "{http://www.w3.org/2000/svg}svg"
    except (OSError, ET.ParseError):
        return False


def render_docs(diagram_dir: str = "docs/diagrams", *, render: bool = True) -> dict:
    source_dir = Path(diagram_dir)
    sources = sorted(source_dir.glob("*.mmd"))
    checked = [
        {"source": source.name, "svg": source.with_suffix(".svg").name,
         "ok": _is_svg(source.with_suffix(".svg"))}
        for source in sources
    ]
    renderer = shutil.which("mmdc")
    if renderer is None:
        local = Path("ui/node_modules/.bin/mmdc")
        if local.is_file():
            renderer = str(local.resolve())
    rendered = []
    if render and renderer:
        out_dir = Path("runs/docs/diagrams")
        out_dir.mkdir(parents=True, exist_ok=True)
        for source in sources:
            target = out_dir / source.with_suffix(".svg").name
            args = [renderer, "-i", str(source), "-o", str(target)]
            config = Path("ui/puppeteer-config.json")
            if config.is_file():
                args.extend(["-p", str(config)])
            try:
                result = subprocess.run(args, capture_output=True, text=True, timeout=40)
                rendered.append({
                    "source": source.name, "ok": result.returncode == 0 and _is_svg(target),
                    "stderr": result.stderr.strip() if result.returncode else "",
                })
            except (OSError, subprocess.TimeoutExpired) as exc:
                rendered.append({"source": source.name, "ok": False, "stderr": str(exc)})
    return {
        "sources": len(sources), "renderer_available": bool(renderer),
        "checked": checked, "rendered": rendered,
        "ok": len(checked) >= 12 and all(r["ok"] for r in [*checked, *rendered]),
    }
