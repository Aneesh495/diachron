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


def build_thesis(thesis_dir: str = "docs/thesis", out_path: str = "docs/thesis/THESIS.md") -> dict:
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

    if parts:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n\n".join(parts))

    # Claim-to-evidence ledger check.
    claims_path = os.path.join(thesis_dir, "claims.json")
    ledger = {"claims": 0, "incomplete": []}
    if os.path.exists(claims_path):
        with open(claims_path, encoding="utf-8") as f:
            claims = json.load(f)
        ledger["claims"] = len(claims)
        for c in claims:
            if not all(k in c and c[k] for k in ("claim", "hypothesis", "evidence", "status")):
                ledger["incomplete"].append(c.get("claim", "<unnamed>"))

    return {
        "gate": "G14",
        "chapters_present": present,
        "chapters_missing": missing,
        "assembled": out_path if parts else None,
        "claim_ledger": ledger,
        "ledger_ok": ledger["claims"] > 0 and not ledger["incomplete"],
    }
