"""Release packaging checks.

`make release` confirms the repository is in a releasable state: the extension
imports, the data manifests exist, the public data release excludes restricted
lexical records, and the working tree is committed. It reports what it finds
rather than performing an irreversible publish; pushing and tagging stay with the
operator.
"""

from __future__ import annotations

import os

from ..util.gitmeta import source_commit, working_tree_clean


def run_release(evidence_dir: str | None = None) -> dict:
    checks = {}
    try:
        from diachron import _core  # noqa: F401

        checks["extension_imports"] = True
    except Exception:  # noqa: BLE001
        checks["extension_imports"] = False

    checks["manifests_present"] = os.path.isdir("data/manifests") and bool(
        [f for f in os.listdir("data/manifests") if f.endswith(".json")]
    )
    # The public release ships only synthetic/derived data and acquisition
    # recipes; confirm no restricted raw lexical dumps are tracked.
    import subprocess

    tracked = subprocess.check_output(["git", "ls-files", "-z"]).decode().split("\0")
    checks["no_restricted_raw_tracked"] = not any(
        name.startswith(("data/raw/", "data/cache/")) or name.endswith((".parquet", ".arrow"))
        for name in tracked
    )
    checks["working_tree_clean"] = working_tree_clean()
    checks["source_commit"] = source_commit()

    checks["source_package_checks_passed"] = all(
        value for value in checks.values() if isinstance(value, bool)
    )
    full_verified = False
    if evidence_dir is not None:
        from pathlib import Path
        from ..verify.verifier import _read_json, verify_evidence

        try:
            record = _read_json(Path(evidence_dir) / "ACCEPTANCE.json", 4 * 1024 * 1024)
            if record["profile"] == "full" and record["full_profile_executed"] and record["ready"]:
                full_verified = verify_evidence(evidence_dir, expected_commit=checks["source_commit"])["ok"]
        except (OSError, KeyError, TypeError, ValueError):
            full_verified = False
    checks["ready"] = full_verified and checks["source_package_checks_passed"]
    if not checks["ready"]:
        checks["release_blocker"] = "A verified full execution is required for release."
    return checks
