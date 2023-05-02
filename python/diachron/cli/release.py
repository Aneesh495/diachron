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


def run_release() -> dict:
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
    restricted = []
    for root, _dirs, files in os.walk("data"):
        for f in files:
            if f.endswith(".zip") or "raw" in root.split(os.sep):
                restricted.append(os.path.join(root, f))
    checks["no_restricted_raw_tracked"] = len(restricted) == 0
    checks["working_tree_clean"] = working_tree_clean()
    checks["source_commit"] = source_commit()

    checks["ready"] = all(
        v is True for k, v in checks.items() if isinstance(v, bool)
    )
    return checks
