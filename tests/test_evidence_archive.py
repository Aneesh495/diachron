"""Only declared regular files with exact identities survive archive restoration."""

import hashlib
import io
import json
import tarfile

import pytest

from diachron.io.evidence_archive import restore_evidence
from diachron.io.types import DataError


def archive(tmp_path, extra=None, altered=False):
    path = tmp_path / 'evidence.tar.gz'
    payload = b'{"posterior":"saved"}'
    inventory = {'bundle.json': {'sha256': hashlib.sha256(payload).hexdigest(), 'bytes': len(payload)}}
    with tarfile.open(path, 'w:gz') as target:
        for name, content in [('case/ARCHIVE_CONTENTS.json', json.dumps(inventory).encode()),
                              ('case/bundle.json', payload[:-1] + b'X' if altered else payload)]:
            member = tarfile.TarInfo(name)
            member.size = len(content)
            target.addfile(member, io.BytesIO(content))
        if extra is not None:
            member = tarfile.TarInfo(extra)
            member.type = tarfile.SYMTYPE if extra == 'case/link' else tarfile.REGTYPE
            member.linkname = '/etc/passwd'
            target.addfile(member)
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def test_verified_restore_and_external_digest(tmp_path):
    path, digest = archive(tmp_path)
    out = tmp_path / 'restored'
    report = restore_evidence(path, out, digest)
    assert report['restored_files'] == 1
    assert json.loads((out / 'bundle.json').read_text()) == {'posterior': 'saved'}
    with pytest.raises(DataError, match='compressed archive digest'):
        restore_evidence(path, tmp_path / 'wrong', '0' * 64)


@pytest.mark.parametrize('extra', ['../escape', '/tmp/escape', 'case/../../escape', 'case/link', 'case/undeclared'])
def test_restore_rejects_escape_links_and_uninventoried_files(tmp_path, extra):
    path, digest = archive(tmp_path, extra)
    with pytest.raises(DataError):
        restore_evidence(path, tmp_path / 'out', digest)
    assert not (tmp_path / 'out').exists()


def test_restore_rejects_changed_member_and_budget(tmp_path):
    path, digest = archive(tmp_path, altered=True)
    with pytest.raises(DataError, match='restored file digest'):
        restore_evidence(path, tmp_path / 'out', digest)
    with pytest.raises(DataError, match='uncompressed budget'):
        restore_evidence(path, tmp_path / 'out', digest, max_bytes=10)
