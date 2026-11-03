#!/usr/bin/env python3
"""Inspect wheel/sdist contents and metadata without executing the artifacts."""
import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import zipfile


def check(path, public=False):
    if path.name.endswith('.whl'):
        with zipfile.ZipFile(path) as archive:
            files = {n: archive.read(n) for n in archive.namelist() if not n.endswith('/')}
        metadata_name = next(n for n in files if n.endswith('.dist-info/METADATA'))
        assert 'rivin_jones/colored.py' in files and 'rivin_jones/r_matrix.py' in files
        assert 'rivin_jones/tensor.py' in files
        assert any(n.startswith('rivin_jones/_tl_tensor.') and n.endswith(('.so', '.pyd')) for n in files)
        assert not any(n.startswith('tl_tensor/') for n in files), 'old import namespace leaked'
        kind = 'wheel'
    else:
        with tarfile.open(path, 'r:gz') as archive:
            members = archive.getmembers()
            modes = {m.name: m.mode for m in members}
            assert not any(m.issym() or m.islnk() for m in members), 'source archive contains links'
            files = {m.name: archive.extractfile(m).read() for m in members if m.isfile()}
        metadata_name = next(n for n in files if n.endswith('/PKG-INFO'))
        prefix = metadata_name.removesuffix('/PKG-INFO')+'/'
        for required in ('Cargo.lock', 'src/lib.rs', 'python/rivin_jones/r_matrix.py',
                         'python/rivin_jones/tensor.py', 'tests/test_tensor.py',
                         'tests/test_invariants.py', 'tests/data/colored.json'):
            assert prefix+required in files, required
        kind = 'sdist'
        if public:
            config = files[prefix+'.cargo/config.toml'].decode()
            assert 'directory = "vendor"' in config, 'vendor configuration is not portable'
            assert prefix+'docs/REBUILDING.md' in files
            crates = 0
            for block in files[prefix+'Cargo.lock'].decode().split('[[package]]')[1:]:
                if 'source = "registry+' not in block:
                    continue
                name = re.search(r'^name = "([^"]+)"', block, re.M).group(1)
                version = re.search(r'^version = "([^"]+)"', block, re.M).group(1)
                checksum = re.search(r'^checksum = "([^"]+)"', block, re.M).group(1)
                base = prefix+'vendor/'+name+'-'+version+'/'
                manifest = json.loads(files[base+'.cargo-checksum.json'])
                assert manifest['package'] == checksum, name
                for source, digest in manifest['files'].items():
                    assert hashlib.sha256(files[base+source]).hexdigest() == digest, base+source
                crates += 1
            assert crates >= 70, 'locked dependency sources missing'
            configure = prefix+'vendor/gmp-mpfr-sys-1.7.1/gmp-6.3.0-c/configure'
            assert configure in files
            assert modes[configure] & 0o111, 'GMP configure lost its executable mode'
    for name in files:
        p = PurePosixPath(name)
        assert not p.is_absolute() and '..' not in p.parts, name
        assert not any(part.startswith('._') or part in ('.git', '__pycache__', '.DS_Store') for part in p.parts), name
        assert not name.endswith(('.pyc', '.prof')), name
    metadata = BytesParser().parsebytes(files[metadata_name])
    assert metadata['Name'] == 'rivin-jones' and metadata['Version'] == '0.1.0'
    assert set(metadata['Requires-Python'].replace(' ', '').split(',')) == {'>=3.10', '<3.14'}
    assert 'tensor' in metadata.get_all('Provides-Extra', [])
    assert any('numpy' in item and 'tensor' in item for item in metadata.get_all('Requires-Dist', []))
    private = 'Private :: Do Not Upload' in metadata.get_all('Classifier', [])
    if 'Redistribution terms pending' in metadata.get('License', ''):
        assert private, 'pending release lost its upload guard'
    if public:
        assert not private, 'candidate is marked private'
        assert metadata.get('License-Expression') == 'GPL-3.0-or-later', 'incorrect project license'
        assert 'LICENSE' in metadata.get_all('License-File', []), 'GPL text is missing'
    assert any('LICENSE_STATUS.md' in n for n in files)
    return dict(path=path.name, kind=kind, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), files=len(files), status='passed')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('artifacts', nargs='+', type=Path)
    p.add_argument('--public', action='store_true', help='also reject private release candidates')
    a = p.parse_args()
    rows = [check(f, a.public) for f in a.artifacts if not f.name.startswith('._')]
    assert rows
    print(json.dumps(rows, indent=2))


if __name__ == '__main__': main()
