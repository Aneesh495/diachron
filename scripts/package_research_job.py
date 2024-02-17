"""Package a completed or failed research job with identity and source rights."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--name', required=True)
parser.add_argument('--directory', type=Path, help='existing run directory, defaults to runs/NAME')
args = parser.parse_args()
if not args.name.replace('_', '').isalnum():
    raise ValueError('invalid job archive name')
root = args.directory or Path('runs') / args.name
if args.directory is not None and not root.is_dir():
    raise FileNotFoundError(root)
root.mkdir(parents=True, exist_ok=True)
rights = root / 'provenance'
rights.mkdir(exist_ok=True)
for source in ('docs/DATA_CARD.md', 'docs/OBSERVED_DATA.md', 'docs/EVIDENCE_ARCHIVES.md', 'LICENSE'):
    shutil.copy2(source, rights / Path(source).name)
shutil.copytree('data/manifests', rights / 'manifests', dirs_exist_ok=True)
inventory = root / 'ARCHIVE_CONTENTS.json'
identities = {}
for path in sorted(root.rglob('*')):
    if path.is_symlink() or not (path.is_file() or path.is_dir()):
        raise ValueError(f'archive inputs must be regular files or directories: {path}')
    if path == inventory or path.is_dir():
        continue
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    identities[path.relative_to(root).as_posix()] = {'sha256': digest, 'bytes': path.stat().st_size}
inventory.write_text(json.dumps(identities, indent=2) + '\n')
archive = Path('runs') / f'{args.name}.tar.gz'
with tarfile.open(archive, 'w:gz', compresslevel=6) as target:
    target.add(root, arcname=args.name)
with archive.open('rb') as stream:
    digest = hashlib.file_digest(stream, 'sha256').hexdigest()
archive.with_suffix('.gz.sha256').write_text(f'{digest}  {archive.name}\n')
print(json.dumps({'archive': str(archive), 'bytes': archive.stat().st_size, 'sha256': digest}))
