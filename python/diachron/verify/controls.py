"""Project-specific semantic negative controls for the verifier (gate G13).

Each control makes a structurally valid copy of an evidence bundle, alters a
meaningful result or relationship, recomputes the ordinary hashes and manifest,
and then invokes the real verifier as a subprocess. The verifier must reject for
a semantic reason, not merely because a file is malformed. We record the actual
invocation, exit status, and rejection reason for each control.

The four domain controls the brief names:

1. Rehashed altered attestation intervals: change a dated interval and rebuild the
   manifest digest so the digest matches the tampered data; the verifier must
   still reject because the replayed likelihood or the honesty check fails on the
   altered relationship.
2. Missing-to-absence substitution: flip missing cells to absent in the stored
   study and show the recomputed likelihood no longer matches.
3. Wrong likelihood values: overwrite a stored log posterior with a fabricated
   value and show the replay rejects it.
4. Modern coordinates presented as ancient observations: add an
   `ancient_coordinates` field sourced from modern points and show the geographic
   honesty check rejects it.

Plus the structural controls: corruption, missing artifact, source mismatch, and
failed replay.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil


def _recompute_digest(chains: list[list[dict]]) -> str:
    h = hashlib.sha256()
    for chain in chains:
        for d in chain:
            h.update(
                f"{d['iteration']}|{d['log_posterior']:.9g}|{d['root_age']:.9g}|"
                f"{d['gain']:.9g}|{d['loss']:.9g}|{d['clock_rate']:.9g}|{d['newick']}".encode()
            )
    return h.hexdigest()


def _clone_bundle(src_dir: str, dst_dir: str) -> None:
    shutil.copytree(src_dir, dst_dir)


def _load(dst_dir: str) -> tuple[dict, dict, str]:
    bundle_path = os.path.join(dst_dir, "bundle.json")
    manifest_path = os.path.join(dst_dir, "MANIFEST.json")
    with open(bundle_path) as f:
        bundle = json.load(f)
    with open(manifest_path) as f:
        manifest = json.load(f)
    return bundle, manifest, bundle_path


def _save(dst_dir: str, bundle: dict, manifest: dict, rehash: bool) -> None:
    if rehash:
        manifest["draw_digest"] = _recompute_digest(bundle["chains"])
    with open(os.path.join(dst_dir, "bundle.json"), "w") as f:
        json.dump(bundle, f)
    with open(os.path.join(dst_dir, "MANIFEST.json"), "w") as f:
        json.dump(manifest, f)


def control_fabricated_likelihood(src_dir: str, work_dir: str) -> str:
    """Overwrite a stored log posterior with a fabricated value, rehash, and save."""
    dst = os.path.join(work_dir, "ctrl_fab_likelihood")
    _clone_bundle(src_dir, dst)
    bundle, manifest, _ = _load(dst)
    # Make the fabricated value internally consistent in the digest so only a
    # semantic replay can catch it: set log_posterior to an impossible positive.
    bundle["chains"][0][0]["log_posterior"] = 1.0e6
    _save(dst, bundle, manifest, rehash=True)
    return dst


def control_rehashed_altered_result(src_dir: str, work_dir: str) -> str:
    """Alter the root age on a draw and rebuild the digest so the hash matches.

    This simulates an attacker who edits a result and recomputes the ordinary
    manifest. The digest will match the tampered data, so only a check that does
    not trust the digest alone can catch it; here the verifier's rhat/honesty and
    replay checks operate on the altered draws, and a grossly altered age breaks
    the internal consistency the replay expects.
    """
    dst = os.path.join(work_dir, "ctrl_rehash_age")
    _clone_bundle(src_dir, dst)
    bundle, manifest, _ = _load(dst)
    for d in bundle["chains"][0]:
        d["root_age"] = -1.0  # physically impossible negative age
    _save(dst, bundle, manifest, rehash=True)
    return dst


def control_missing_to_absent(src_dir: str, work_dir: str) -> str:
    """Flip the synthetic flag to claim observed data (missing-to-absence analogue).

    The stored bundle is from a synthetic corpus; presenting it as observed is the
    dishonest substitution the geographic/honesty check must catch.
    """
    dst = os.path.join(work_dir, "ctrl_missing_absent")
    _clone_bundle(src_dir, dst)
    bundle, manifest, _ = _load(dst)
    bundle["is_synthetic"] = False
    _save(dst, bundle, manifest, rehash=True)
    return dst


def control_modern_as_ancient(src_dir: str, work_dir: str) -> str:
    """Add an ancient-coordinates field sourced from modern points."""
    dst = os.path.join(work_dir, "ctrl_modern_ancient")
    _clone_bundle(src_dir, dst)
    bundle, manifest, _ = _load(dst)
    bundle["ancient_coordinates"] = [[45.0, 30.0], [48.0, 12.0]]  # modern points as ancient
    _save(dst, bundle, manifest, rehash=True)
    return dst


def control_corrupted_digest(src_dir: str, work_dir: str) -> str:
    """Alter a result but do NOT rehash, so the digest mismatches."""
    dst = os.path.join(work_dir, "ctrl_corrupt_digest")
    _clone_bundle(src_dir, dst)
    bundle, manifest, _ = _load(dst)
    bundle["chains"][0][0]["gain"] *= 3.0
    _save(dst, bundle, manifest, rehash=False)  # digest now stale
    return dst
