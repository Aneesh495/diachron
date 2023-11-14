"""Diachron research CLI.

A single entry point dispatches subcommands for the whole research workflow:
inspect the host (`doctor`), build data (`data`), run the likelihood oracle
(`likelihood-check`), sample (`sample`), diagnose (`diagnose`), run posterior
predictive and sensitivity checks, run the recovery study, build the thesis, and
drive the end-to-end demonstration. Long campaigns are bounded by explicit
iteration and chain arguments; nothing here spawns unbounded work.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import json
import os
import platform
import shutil
import sys
import time

from .. import __version__


def _print(obj) -> None:
    if isinstance(obj, (dict, list)):
        print(json.dumps(obj, indent=2))
    else:
        print(obj)


def cmd_doctor(args) -> int:
    """Report toolchain, CPU, storage, and data availability (D30)."""
    info = {
        "diachron_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
    }
    try:
        import numpy

        info["numpy"] = numpy.__version__
    except ImportError:
        info["numpy"] = "MISSING"
    try:
        from diachron import _core

        info["core_extension"] = "ok"
        info["core_haversine_check"] = round(_core.haversine_km(0, 0, 0, 1), 3)
    except Exception as exc:  # noqa: BLE001
        info["core_extension"] = f"MISSING: {exc}"
    usage = shutil.disk_usage(".")
    info["disk_free_gb"] = round(usage.free / 1e9, 1)
    info["data_manifests"] = sorted(
        f for f in os.listdir("data/manifests") if f.endswith(".json")
    ) if os.path.isdir("data/manifests") else []
    _print(info)
    return 0


def cmd_data(args) -> int:
    """Build the reference corpora and write manifests."""
    from ..io.corpus import build_primary_corpus, build_second_corpus, write_manifest

    os.makedirs("data/manifests", exist_ok=True)
    prim = build_primary_corpus()
    m1 = write_manifest(prim, "data/manifests/primary.json")
    sec = build_second_corpus()
    m2 = write_manifest(sec, "data/manifests/second.json")
    _print({"primary": m1["dimensions"], "second": m2["dimensions"]})
    return 0


def cmd_cldf(args) -> int:
    """Ingest a CLDF wordlist directory into typed tables."""
    from ..io.cldf import load_cldf_wordlist
    from ..io.types import DataError

    if args.action == "freeze":
        from pathlib import Path
        from ..io.freeze import load_study, write_study

        if not all((args.study, args.source, args.catalog, args.out)):
            print("provide --study, --source, --catalog and --out", file=sys.stderr)
            return 2
        try:
            dataset, crosswalk, report = load_study(
                args.study, Path(args.source), Path(args.catalog),
            )
            write_study(Path(args.out), dataset, crosswalk, report)
        except (DataError, OSError, ValueError) as exc:
            print(str(exc), file=sys.stderr)
            return 2
        _print(report)
        return 0

    if not args.metadata:
        print("provide --metadata path to a CLDF *-metadata.json", file=sys.stderr)
        return 2
    ds = load_cldf_wordlist(args.metadata)
    _print(
        {
            "name": ds.name,
            "license": ds.license,
            "languages": ds.n_eligible_languages,
            "concepts": ds.n_eligible_concepts,
            "lexemes": len(ds.lexemes),
            "judgments": len(ds.judgments),
        }
    )
    return 0


def cmd_attestations(args) -> int:
    """Summarize the attestation corpus (observation types, dating spread)."""
    if args.action == "edh":
        from pathlib import Path
        from ..io.epidoc import load_edh
        from ..io.types import DataError

        if not all((args.source, args.geography, args.out)):
            print("provide --source, --geography and --out", file=sys.stderr)
            return 2
        try:
            report = load_edh(Path(args.source), Path(args.geography))
            out = Path(args.out)
            out.parent.mkdir(parents=True, exist_ok=True)
            with out.open("x", encoding="utf-8") as file:
                json.dump(report, file, indent=2, ensure_ascii=False, allow_nan=False)
        except (DataError, OSError, ValueError) as exc:
            print(str(exc), file=sys.stderr)
            return 2
        _print({key: value for key, value in report.items() if key != "attestations"})
        return 0

    from ..io.corpus import build_primary_corpus

    corpus = build_primary_corpus()
    types: dict[str, int] = {}
    for a in corpus.attestations:
        types[a.observation_type] = types.get(a.observation_type, 0) + 1
    _print(
        {
            "n_attestations": len(corpus.attestations),
            "observation_types": types,
            "example": {
                "id": corpus.attestations[0].id,
                "age_low": corpus.attestations[0].age_low,
                "age_high": corpus.attestations[0].age_high,
                "source_ref": corpus.attestations[0].source_ref,
            },
        }
    )
    return 0


def cmd_likelihood_check(args) -> int:
    """Run the pruning-vs-enumeration oracle (gate G02)."""
    from .checks import run_likelihood_oracle

    result = run_likelihood_oracle(n_cases=args.cases, seed=args.seed)
    _print(result)
    return 0 if result["passed"] else 1


def cmd_sample(args) -> int:
    """Run a bounded MCMC campaign and save the posterior bundle."""
    from ..io.types import DataError
    from .sampling import run_sample

    try:
        result = run_sample(
            model_name=args.model,
            n_chains=args.chains,
            warmup=args.warmup,
            iterations=args.iterations,
            thin=args.thin,
            out_dir=args.out,
            profile=args.profile,
            resume=args.resume,
            jobs=args.jobs,
            geography_error=args.geography_error,
            observed_study=args.study,
            observed_source=args.source,
            catalog=args.catalog,
            edh_source=args.edh_source,
            edh_geography=args.edh_geography,
        )
    except (DataError, OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    _print(result)
    return 0


def cmd_diagnose(args) -> int:
    """Compute convergence and mixing diagnostics for a saved bundle."""
    from .sampling import run_diagnose

    result = run_diagnose(args.bundle)
    _print(result)
    return 0


def cmd_ppc(args) -> int:
    from ..ppc.run import run_posterior_predictive

    result = run_posterior_predictive(args.bundle)
    _print(result)
    return 0


def cmd_sensitivity(args) -> int:
    from ..sensitivity.run import run_sensitivity

    result = run_sensitivity(out_dir=args.out)
    _print(result)
    return 0


def cmd_recovery(args) -> int:
    from ..sensitivity.recovery import run_recovery_study

    result = run_recovery_study(n_datasets=args.datasets, out_dir=args.out)
    _print(result)
    return 1 if result["failed"] else 0


def cmd_predict(args) -> int:
    from ..ppc.heldout import run_heldout_prediction

    result = run_heldout_prediction(out_dir=args.out)
    _print(result)
    return 0


def cmd_thesis(args) -> int:
    from ..cli.thesis import build_thesis

    result = build_thesis()
    _print(result)
    return 0 if result["ledger_ok"] and result["rendered_html"] else 1


def cmd_demo(args) -> int:
    from .demo import run_demo

    return run_demo(out_dir=args.out)


def cmd_bench(args) -> int:
    from .checks import run_bench

    _print(run_bench())
    return 0


def cmd_stress(args) -> int:
    from .checks import run_stress

    _print(run_stress())
    return 0


def cmd_docs(args) -> int:
    from .docs import render_docs

    result = render_docs(render=args.action == "render")
    _print(result)
    return 0 if result["ok"] else 1


def cmd_acceptance(args) -> int:
    from .acceptance import run_acceptance

    return run_acceptance(run_id=args.run_id, profile=args.profile)


def cmd_verify(args) -> int:
    from ..verify.verifier import verify_evidence

    result = verify_evidence(
        args.evidence, expected_commit=args.expected_commit, observed_source=args.source,
        catalog=args.catalog, edh_source=args.edh_source, edh_geography=args.edh_geography,
    )
    _print(result)
    return 0 if result["ok"] else 1


def cmd_release(args) -> int:
    from .release import run_release

    _print(run_release())
    return 0


def cmd_export_workbench(args) -> int:
    from .workbench import export_workbench

    bundle_dirs = {}
    for spec in args.bundle:
        name, _, path = spec.partition("=")
        bundle_dirs[name] = path
    result = export_workbench(bundle_dirs, args.out, source=args.source)
    _print(result)
    return 0


def cmd_inspect_attestation(args) -> int:
    from dataclasses import asdict

    from ..infer.artifacts import PosteriorBundle
    from ..io.types import Attestation
    from ..registry import AttestationInspector

    bundle = PosteriorBundle.load(args.bundle)
    records = [Attestation(**raw) for raw in bundle.context.get("attestations", ())]
    try:
        report = AttestationInspector(records).inspect(args.id)
    except KeyError:
        print(f"attestation is absent from this bundle: {args.id}", file=sys.stderr)
        return 2
    _print(asdict(report))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="diachron", description="Diachron research CLI")
    p.add_argument("--version", action="version", version=f"diachron {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("doctor").set_defaults(func=cmd_doctor)

    dp = sub.add_parser("data")
    dp.add_argument("action", nargs="?", default="build")
    dp.set_defaults(func=cmd_data)

    cp = sub.add_parser("cldf")
    cp.add_argument("action", nargs="?", choices=("ingest", "freeze"), default="ingest")
    cp.add_argument("--metadata", default=None)
    cp.add_argument("--study", choices=("iecor", "abvd"))
    cp.add_argument("--source", help="root of the pinned source checkout")
    cp.add_argument("--catalog", help="pinned Glottolog languoid ZIP")
    cp.add_argument("--out", help="new directory for the frozen Arrow tables")
    cp.set_defaults(func=cmd_cldf)

    ap = sub.add_parser("attestations")
    ap.add_argument("action", nargs="?", choices=("build", "edh"), default="build")
    ap.add_argument("--source", help="directory with the pinned EDH XML files")
    ap.add_argument("--geography", help="pinned EDH GeoJSON gazetteer")
    ap.add_argument("--out", help="new attestation JSON file")
    ap.set_defaults(func=cmd_attestations)

    lp = sub.add_parser("likelihood-check")
    lp.add_argument("--cases", type=int, default=1000)
    lp.add_argument("--seed", type=int, default=0)
    lp.set_defaults(func=cmd_likelihood_check)

    sp = sub.add_parser("sample")
    sp.add_argument("--model", default="dated")
    sp.add_argument("--chains", type=int, default=4)
    sp.add_argument("--warmup", type=int, default=1500)
    sp.add_argument("--iterations", type=int, default=6000)
    sp.add_argument("--thin", type=int, default=10)
    sp.add_argument("--out", default="runs/sample")
    sp.add_argument("--profile", choices=("cpu", "full"), default="cpu")
    sp.add_argument("--resume", action="store_true")
    sp.add_argument("--jobs", type=int, default=1)
    sp.add_argument("--geography-error", type=float, help="fixed planar reference error in degrees")
    sp.add_argument("--study", choices=("iecor", "abvd"))
    sp.add_argument("--source", type=Path, help="root of pinned CLDF source checkout")
    sp.add_argument("--catalog", type=Path, help="pinned Glottolog archive")
    sp.add_argument("--edh-source", type=Path, help="directory containing pinned EDH XML files")
    sp.add_argument("--edh-geography", type=Path, help="pinned EDH gazetteer")
    sp.set_defaults(func=cmd_sample)

    gp = sub.add_parser("diagnose")
    gp.add_argument("--bundle", default="runs/sample/dated")
    gp.set_defaults(func=cmd_diagnose)

    pp = sub.add_parser("ppc")
    pp.add_argument("--bundle", default="runs/sample/dated")
    pp.set_defaults(func=cmd_ppc)

    sn = sub.add_parser("sensitivity")
    sn.add_argument("--out", default="runs/sensitivity")
    sn.set_defaults(func=cmd_sensitivity)

    rp = sub.add_parser("recovery")
    rp.add_argument("--datasets", type=int, default=200)
    rp.add_argument("--out", default="runs/recovery")
    rp.set_defaults(func=cmd_recovery)

    pr = sub.add_parser("predict")
    pr.add_argument("--out", default="runs/predict")
    pr.set_defaults(func=cmd_predict)

    tp = sub.add_parser("thesis")
    tp.add_argument("action", nargs="?", default="build")
    tp.set_defaults(func=cmd_thesis)

    mp = sub.add_parser("demo")
    mp.add_argument("--out", default="runs/demo")
    mp.set_defaults(func=cmd_demo)

    sub.add_parser("bench").set_defaults(func=cmd_bench)
    sub.add_parser("stress").set_defaults(func=cmd_stress)

    dc = sub.add_parser("docs")
    dc.add_argument("action", nargs="?", choices=("render", "check"), default="render")
    dc.set_defaults(func=cmd_docs)

    acp = sub.add_parser("acceptance")
    acp.add_argument("--run-id", dest="run_id", required=True)
    acp.add_argument("--profile", choices=("cpu", "full"), default="full")
    acp.set_defaults(func=cmd_acceptance)

    vp = sub.add_parser("verify")
    vp.add_argument("--evidence", required=True)
    vp.add_argument("--expected-commit", help="trusted executable source commit")
    vp.add_argument("--source", type=Path, help="trusted pinned lexical inputs")
    vp.add_argument("--catalog", type=Path, help="trusted Glottolog archive")
    vp.add_argument("--edh-source", type=Path, help="trusted EDH XML inputs")
    vp.add_argument("--edh-geography", type=Path, help="trusted EDH gazetteer")
    vp.set_defaults(func=cmd_verify)

    sub.add_parser("release").set_defaults(func=cmd_release)

    isp = sub.add_parser("inspect-attestation")
    isp.add_argument("--bundle", required=True)
    isp.add_argument("--id", required=True)
    isp.set_defaults(func=cmd_inspect_attestation)

    ewp = sub.add_parser("export-workbench")
    ewp.add_argument("--bundle", action="append", default=[],
                     help="model=path, repeatable (e.g. dated=runs/demo/dated)")
    ewp.add_argument("--out", default="ui/public/workbench.json")
    ewp.add_argument("--source", type=Path, help="optional pinned wordlist for archived display labels")
    ewp.set_defaults(func=cmd_export_workbench)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    start = time.time()
    rc = args.func(args)
    if os.environ.get("DIACHRON_TIMING"):
        print(f"[{args.command} took {time.time() - start:.1f}s]", file=sys.stderr)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
