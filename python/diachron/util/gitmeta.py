"""Record the executable source commit so evidence points at exact code.

Every posterior bundle and evidence directory stores the commit it was produced
from. The verifier later checks that commit against the trusted checkout rather
than trusting the bundle's own claim, so this function only reports what it can
read from the working tree; it never asserts the commit is published.
"""

from __future__ import annotations

import subprocess


def source_commit() -> str:
    """Return the current git HEAD short hash, or 'unknown' outside a repo."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def working_tree_clean() -> bool:
    """True if there are no uncommitted changes to tracked files."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        return out.returncode == 0 and out.stdout.strip() == ""
    except (OSError, subprocess.SubprocessError):
        return False
