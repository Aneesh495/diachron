"""Restore bounded evidence archives without trusting tar paths or links."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile

from .types import DataError

_MAX_BYTES = 2_000_000_000
_MAX_MEMBERS = 20000
_MAX_FILE = 128 * 1024 * 1024


def _relative(name):
    path = PurePosixPath(name)
    if path.is_absolute() or not path.parts or any(part in {'.', '..'} for part in path.parts) or '\\' in name:
        raise DataError('unsafe path in evidence archive')
    return path


def restore_evidence(archive: Path, out: Path, expected_sha256: str, max_bytes=_MAX_BYTES):
    if not 0 < max_bytes <= _MAX_BYTES or not archive.is_file() or archive.stat().st_size > _MAX_BYTES:
        raise DataError('invalid archive size or restore budget')
    if len(expected_sha256) != 64 or any(c not in '0123456789abcdef' for c in expected_sha256):
        raise DataError('provide a trusted lowercase SHA-256 archive digest')
    with archive.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    if digest != expected_sha256:
        raise DataError('compressed archive digest differs')
    if out.exists():
        raise FileExistsError('evidence restore output exists')
    inventory = None
    inventory_bytes = None
    top = None
    member_names = set()
    total = 0
    with tarfile.open(archive, 'r|gz') as source:
        for index, member in enumerate(source):
            if index >= _MAX_MEMBERS:
                raise DataError('too many archive members')
            relative = _relative(member.name)
            if member.name in member_names:
                raise DataError('duplicate archive member')
            member_names.add(member.name)
            if top is None:
                top = relative.parts[0]
            if relative.parts[0] != top or not (member.isdir() or member.isfile()):
                raise DataError('archive contains links, special files or multiple roots')
            if member.isfile():
                total += member.size
                if member.size < 0 or member.size > _MAX_FILE or total > max_bytes:
                    raise DataError('archive exceeds its uncompressed budget')
                if relative == PurePosixPath(top) / 'ARCHIVE_CONTENTS.json':
                    if member.size > 8 * 1024 * 1024:
                        raise DataError('oversized archive inventory')
                    inventory_bytes = source.extractfile(member).read(member.size + 1)
                    inventory = json.loads(inventory_bytes)
    if not isinstance(inventory, dict) or not inventory:
        raise DataError('archive inventory missing or empty')
    expected_files = set()
    for name, identity in inventory.items():
        relative = _relative(name)
        if relative.parts[0] == top or name == 'ARCHIVE_CONTENTS.json':
            raise DataError('inventory paths must be relative to the archive root')
        if type(identity['bytes']) is not int or not 0 <= identity['bytes'] <= _MAX_FILE:
            raise DataError('invalid inventory byte count')
        expected_files.add(str(PurePosixPath(top) / relative))
    actual_files = set()
    out.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(out.parent).free < total + 50_000_000:
        raise DataError('insufficient space for bounded evidence restoration')
    with tempfile.TemporaryDirectory(dir=out.parent, prefix='.restore-') as temporary:
        stage = Path(temporary) / 'evidence'
        stage.mkdir()
        with tarfile.open(archive, 'r|gz') as source:
            for member in source:
                relative = _relative(member.name)
                target = stage.joinpath(*relative.parts[1:])
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                if relative == PurePosixPath(top) / 'ARCHIVE_CONTENTS.json':
                    continue
                if member.name not in expected_files:
                    raise DataError('archive contains an undeclared file')
                identity = inventory[str(PurePosixPath(*relative.parts[1:]))]
                if member.size != identity['bytes']:
                    raise DataError('restored byte count differs from inventory')
                target.parent.mkdir(parents=True, exist_ok=True)
                hashed = hashlib.sha256()
                with source.extractfile(member) as incoming, target.open('xb') as output:
                    while chunk := incoming.read(1024 * 1024):
                        hashed.update(chunk)
                        output.write(chunk)
                if hashed.hexdigest() != identity['sha256']:
                    raise DataError('restored file digest differs from inventory')
                actual_files.add(member.name)
        if actual_files != expected_files:
            raise DataError('archive omits declared inventory files')
        (stage / 'ARCHIVE_CONTENTS.json').write_bytes(inventory_bytes)
        stage.rename(out)
    return {'archive_sha256': digest, 'restored_files': len(actual_files),
            'uncompressed_bytes': total, 'out': str(out)}
