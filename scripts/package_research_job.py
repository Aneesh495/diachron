"""Package a completed or failed research job with identity and source rights."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--name', required=True)
args = parser.parse_args()
if not args.name.replace('_', '').isalnum():
    raise ValueError('invalid job archive name')
root = Path('runs') / args.name
root.mkdir(parents=True, exist_ok=True)
rights = root / 'provenance'
rights.mkdir(exist_ok=True)
for source in ('docs/DATA_CARD.md', 'docs/OBSERVED_DATA.md', 'docs/EVIDENCE_ARCHIVES.md', 'LICENSE'):
    shutil.copy2(source, rights / Path(source).name)
shutil.copytree('data/manifests', rights / 'manifests', dirs_exist_ok=True)
identities = {str(path.relative_to(root)): {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'bytes': path.stat().st_size}
              for path in sorted(root.rglob('*')) if path.is_file() and path.name != 'ARCHIVE_CONTENTS.json'}
(root / 'ARCHIVE_CONTENTS.json').write_text(json.dumps(identities, indent=2) + '\n')
archive = Path('runs') / f'{args.name}.tar.gz'
with tarfile.open(archive, 'w:gz', compresslevel=6) as target:
    target.add(root, arcname=args.name)
digest = hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_suffix('.gz.sha256').write_text(f'{digest}  {archive.name}\n')
print(json.dumps({'archive': str(archive), 'bytes': archive.stat().st_size, 'sha256': digest}))
