"""Acquisition rejects oversized, redirected and altered source bytes."""

import hashlib
import io

import pytest

from diachron.io.acquire import _download
from diachron.io.types import DataError


class Response(io.BytesIO):
    def geturl(self):
        return 'https://raw.githubusercontent.com/example/data'


def test_download_checks_bytes_and_refuses_existing_corruption(tmp_path, monkeypatch):
    import diachron.io.acquire as acquire

    payload = b'fixed source bytes'
    expected = {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
    monkeypatch.setattr(acquire, 'urlopen', lambda *a, **k: Response(payload))
    item = ('cldf/table.csv', 'https://raw.githubusercontent.com/example/data', expected)
    assert _download(tmp_path, item) == 'cldf/table.csv'
    assert _download(tmp_path, item) == 'cldf/table.csv'
    (tmp_path / 'cldf/table.csv').write_bytes(payload + b'x')
    with pytest.raises(DataError, match='existing acquisition'):
        _download(tmp_path, item)
    monkeypatch.setattr(acquire, 'urlopen', lambda *a, **k: Response(payload + b'x'))
    with pytest.raises(DataError, match='pinned byte identity'):
        _download(tmp_path, ('other.csv', item[1], expected))
    assert not (tmp_path / 'other.csv').exists()
    assert not list(tmp_path.glob('.download-*'))
    with pytest.raises(DataError, match='escapes'):
        _download(tmp_path, ('../outside.csv', item[1], expected))
    with pytest.raises(DataError, match='outside the declared'):
        _download(tmp_path, ('other.csv', 'https://example.com/data', expected))


def test_declared_git_checkout_conversion_still_requires_exact_hash(tmp_path, monkeypatch):
    import diachron.io.acquire as acquire

    blob = b'id,value\n1,example\n'
    checkout = blob.replace(b'\n', b'\r\n')
    expected = {'bytes': len(checkout), 'sha256': hashlib.sha256(checkout).hexdigest(), 'checkout_eol': 'crlf'}
    monkeypatch.setattr(acquire, 'urlopen', lambda *a, **k: Response(blob))
    _download(tmp_path, ('table.csv', 'https://raw.githubusercontent.com/example/data', expected))
    assert (tmp_path / 'table.csv').read_bytes() == checkout
    monkeypatch.setattr(acquire, 'urlopen', lambda *a, **k: Response(blob.replace(b'example', b'changed')))
    with pytest.raises(DataError, match='pinned byte identity'):
        _download(tmp_path, ('changed.csv', 'https://raw.githubusercontent.com/example/data', expected))
