#!/usr/bin/env python3
"""Verify every public PyPI release file against the tested artifact manifest."""
import hashlib
import json
from pathlib import Path
import sys
import time
import urllib.error
import urllib.request

manifest = json.loads(Path(sys.argv[1]).read_text())
url = 'https://pypi.org/pypi/rivin-jones/' + manifest['version'] + '/json'
for attempt in range(12):
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            release = json.load(response)
        if len(release['urls']) == len(manifest['files']):
            break
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise
    time.sleep(10)
else:
    raise RuntimeError('Public index did not expose the complete release')
remote = {f['filename']: f for f in release['urls']}
assert set(remote) == set(manifest['files']), 'Public file set differs from tested set'
for name, digest in manifest['files'].items():
    item = remote[name]
    assert not item['yanked']
    assert item['digests']['sha256'] == digest, name
    h = hashlib.sha256()
    with urllib.request.urlopen(item['url'], timeout=60) as response:
        for chunk in iter(lambda: response.read(1024*1024), b''):
            h.update(chunk)
    assert h.hexdigest() == digest, name
    print('Verified public download:', name)
Path('public-verification.json').write_text(json.dumps({
    'status': 'passed', 'index': url, 'manifest': manifest,
    'verification': 'All public files downloaded and SHA-256 checked against tested artifacts'
}, indent=2)+'\n')
