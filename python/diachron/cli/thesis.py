"""Assemble the thesis and check its claim-to-evidence ledger (gate G14).

The thesis chapters live as Markdown under `docs/thesis/`. This builder
concatenates them in order into a single document, and validates the argument
ledger: every entry in `docs/thesis/claims.json` must name a claim, the
hypothesis it bears on, the evidence artifact, and a status. The builder reports
any claim missing its evidence pointer so an empirical conclusion can never be
presented without a trace to the evidence behind it.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit, urlunsplit


CHAPTER_ORDER = [
    "00_abstract.md",
    "01_research_questions.md",
    "02_related_work.md",
    "03_evidence_and_annotation.md",
    "04_generative_models.md",
    "05_likelihood_derivations.md",
    "06_sampler_correctness.md",
    "07_experimental_protocol.md",
    "08_synthetic_recovery.md",
    "09_real_data_results.md",
    "10_borrowing.md",
    "11_geography.md",
    "12_sensitivity.md",
    "13_counterarguments.md",
    "14_limitations.md",
]


def build_thesis(thesis_dir: str = "docs/thesis", out_path: str = "runs/thesis/THESIS.md",
                 reference_base: str | None = None) -> dict:
    import markdown
    from markdown.extensions import Extension
    from markdown.treeprocessors import Treeprocessor

    if reference_base:
        base = urlsplit(reference_base)
        if base.scheme != "https" or not base.netloc or base.query or base.fragment:
            raise ValueError("the thesis reference base must be an HTTPS directory URL")

    source_dir = Path(thesis_dir).resolve()
    output_dir = Path(out_path).resolve().parent
    checkout = Path.cwd().resolve()

    class References(Treeprocessor):
        def run(self, root):
            for element in root.iter():
                attribute = {"a": "href", "img": "src"}.get(element.tag)
                if attribute is None or not (value := element.get(attribute)):
                    continue
                url = urlsplit(value)
                if url.scheme or url.netloc or not url.path or url.path.startswith("/"):
                    continue
                source = (source_dir / unquote(url.path)).resolve()
                if reference_base:
                    path = f"{base.path.rstrip('/')}/{quote(source.relative_to(checkout).as_posix())}"
                    rewritten = urlunsplit((base.scheme, base.netloc, path, url.query, url.fragment))
                else:
                    path = quote(Path(os.path.relpath(source, output_dir)).as_posix())
                    rewritten = urlunsplit(("", "", path, url.query, url.fragment))
                element.set(attribute, rewritten)

    class ThesisReferences(Extension):
        def extendMarkdown(self, md):
            md.treeprocessors.register(References(md), "diachron_references", 3)

    rendered = None
    present = []
    missing = []
    parts = []
    for ch in CHAPTER_ORDER:
        path = os.path.join(thesis_dir, ch)
        if os.path.exists(path):
            present.append(ch)
            with open(path, encoding="utf-8") as f:
                parts.append(f.read().rstrip() + "\n")
        else:
            missing.append(ch)

    if parts and not missing:
        combined = "\n\n".join(parts)
        target = Path(out_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(combined, encoding="utf-8")
        body = markdown.markdown(combined, extensions=("fenced_code", "tables", "toc", ThesisReferences()))
        rendered = target.with_suffix(".html")
        rendered.write_text(
            """<!doctype html>
<html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Diachron research report</title>
<style>
body { font: 16px/1.65 Georgia, serif; color: #1c2733; margin: 2rem auto; max-width: 74ch; padding: 0 1rem; }
h1, h2, h3 { line-height: 1.2; }
pre { overflow: auto; background: #eff2f5; padding: 0.8rem; }
td, th { padding: 0.3rem 0.5rem; border-bottom: 1px solid #b9c3ce; }
@media print { h1 { break-before: page; } h1:first-child { break-before: auto; } }
</style>
<main>""" + body + "</main></html>\n",
            encoding="utf-8",
        )

    # Claim-to-evidence ledger check.
    claims_path = os.path.join(thesis_dir, "claims.json")
    ledger = {"claims": 0, "incomplete": [], "unresolved_evidence": []}
    if os.path.exists(claims_path):
        with open(claims_path, encoding="utf-8") as f:
            claims = json.load(f)
        ledger["claims"] = len(claims)
        for c in claims:
            fields = ("claim", "hypothesis", "evidence", "counterevidence", "status", "scope")
            if not all(isinstance(c.get(k), str) and c[k].strip() for k in fields):
                ledger["incomplete"].append(c.get("claim", "<unnamed>"))
                continue
            evidence = Path(c["evidence"].split(":", 1)[0])
            if evidence.is_absolute() or ".." in evidence.parts or not evidence.is_file():
                ledger["unresolved_evidence"].append(c["evidence"])

    return {
        "gate": "G14",
        "chapters_present": present,
        "chapters_missing": missing,
        "assembled": out_path if rendered else None,
        "rendered_html": str(rendered) if rendered else None,
        "claim_ledger": ledger,
        "ledger_ok": bool(ledger["claims"]) and not (
            ledger["incomplete"] or ledger["unresolved_evidence"]
        ),
    }
