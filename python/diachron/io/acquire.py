"""Acquire only the immutable bytes declared in the observed manifests."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import tempfile
import ssl
import time
from urllib.error import URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

from .types import DataError

_MAX_INPUT_BYTES = 250_000_000
_ALLOWED_HOSTS = {"raw.githubusercontent.com", "cdstar.eva.mpg.de"}


def _manifest(name):
    return json.loads(Path(f"data/manifests/{name}.json").read_text(encoding="utf-8"))


def _valid_file(path, expected):
    if not path.is_file() or path.stat().st_size != expected["bytes"]:
        return False
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"]


def _download(root, item):
    relative, url, expected = item
    path = (root / relative).resolve()
    if not path.is_relative_to(root) or path == root:
        raise DataError("acquisition path escapes its output directory")
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in _ALLOWED_HOSTS:
        raise DataError("acquisition URL is outside the declared source hosts")
    if path.exists():
        if not _valid_file(path, expected):
            raise DataError(f"existing acquisition bytes differ: {relative}")
        return relative
    if not 0 < expected["bytes"] <= 64_000_000:
        raise DataError(f"invalid input size: {relative}")
    path.parent.mkdir(parents=True, exist_ok=True)
    context = ssl.create_default_context()
    if not context.cert_store_stats()["x509_ca"]:
        for bundle in (Path("/etc/ssl/cert.pem"), Path("/opt/homebrew/etc/openssl@3/cert.pem")):
            if bundle.is_file():
                context.load_verify_locations(cafile=str(bundle))
                break
    for attempt in range(3):
        fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=".download-")
        try:
            request = Request(url, headers={"User-Agent": "Diachron-reproducible-acquisition/0.1"})
            with urlopen(request, timeout=30, context=context) as response, os.fdopen(fd, "wb") as target:
                fd = None
                final = urlparse(response.geturl())
                if final.scheme != "https" or final.hostname not in _ALLOWED_HOSTS:
                    raise DataError("source redirected outside its allowed host")
                remaining = expected["bytes"] + 1
                while remaining:
                    chunk = response.read(min(1024 * 1024, remaining))
                    if not chunk:
                        break
                    target.write(chunk)
                    remaining -= len(chunk)
                target.flush()
                os.fsync(target.fileno())
            if not _valid_file(Path(temporary), expected) and expected.get("checkout_eol") == "crlf":
                # Reproduce the pinned repository's declared Git checkout rule.
                # The pre-existing observed-byte hash must still match exactly.
                raw = Path(temporary).read_bytes()
                Path(temporary).write_bytes(raw.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
            if not _valid_file(Path(temporary), expected):
                raise DataError(f"download differs from pinned byte identity: {relative}")
            os.replace(temporary, path)
            return relative
        except (URLError, TimeoutError, ConnectionError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
        finally:
            if fd is not None:
                os.close(fd)
            Path(temporary).unlink(missing_ok=True)


def acquire_inputs(out: Path, studies=("iecor",), include_edh=True, workers=4):
    if not 1 <= workers <= 4 or not studies or any(name not in {"iecor", "abvd"} for name in studies):
        raise ValueError("acquisition requires named studies and one to four workers")
    out = out.resolve()
    items = []
    out.mkdir(parents=True, exist_ok=True)
    transport = _manifest("transport")["git_checkout_attributes"]
    for name in sorted(set(studies)):
        manifest = _manifest(name)
        base = f"https://raw.githubusercontent.com/lexibank/{name}/{manifest['git_commit']}/"
        attributes = transport[name]
        attribute_item = (f"{name}/{attributes['path']}", base + attributes["path"],
                          {key: attributes[key] for key in ("bytes", "sha256")})
        _download(out, attribute_item)
        if (out / attribute_item[0]).read_text() != attributes["rule"]:
            raise DataError("source Git checkout attributes differ from the declared transport rule")
        items.append(attribute_item)
        items.extend((f"{name}/{relative}", base + quote(relative), {
            **expected, **({"checkout_eol": "crlf"} if relative.startswith("cldf/") and relative.endswith(".csv") else {}),
        }) for relative, expected in manifest["files"].items())
    catalog = _manifest("glottolog")
    items.append(("glottolog53.zip", catalog["archive_url"], catalog["archive"]))
    if include_edh:
        edh = _manifest("edh")
        base = f"https://raw.githubusercontent.com/epigraphic-database-heidelberg/data/{edh['git_commit']}/"
        items.extend((f"edh/{name}", base + f"inscriptions/1/1/{quote(name)}", expected)
                     for name, expected in edh["files"].items())
        items.append(("edh-geography.json", base + "geography/edhGeographicData.json", edh["geography"]))
    total = sum(expected["bytes"] for _, _, expected in items)
    if total > _MAX_INPUT_BYTES or len({relative for relative, _, _ in items}) != len(items):
        raise DataError("acquisition plan exceeds its byte budget or repeats paths")
    out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        completed = list(pool.map(lambda item: _download(out, item), items))
    report = {
        "schema": 1, "studies": sorted(set(studies)), "edh": include_edh,
        "source_files": len(completed), "total_source_bytes": total,
        "maximum_source_bytes": _MAX_INPUT_BYTES,
        "files": {relative: {"url": url, **expected} for relative, url, expected in items},
    }
    temporary = out / "ACQUISITION.json.tmp"
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(out / "ACQUISITION.json")
    return report
