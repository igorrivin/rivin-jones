#!/usr/bin/env python3
"""Record or verify the exact distribution set produced by successful CI."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['create', 'verify'])
    p.add_argument('directory', type=Path)
    p.add_argument('--commit', required=True)
    p.add_argument('--run-id', required=True)
    a = p.parse_args()
    assert re.fullmatch(r'[0-9a-f]{40}', a.commit)
    assert a.run_id.isdigit()
    files = sorted(x for x in a.directory.iterdir()
                   if x.name.endswith(('.whl', '.tar.gz')) and not x.name.startswith('._'))
    assert len([x for x in files if x.suffix == '.whl']) == 12
    assert len([x for x in files if x.name.endswith('.tar.gz')]) == 1
    manifest = dict(version='0.1.0', source_commit=a.commit, validation_run_id=a.run_id,
                    files={x.name: hashlib.sha256(x.read_bytes()).hexdigest() for x in files})
    target = a.directory / 'SHA256SUMS.json'
    if a.command == 'create':
        target.write_text(json.dumps(manifest, indent=2)+'\n')
    else:
        assert json.loads(target.read_text()) == manifest, 'Artifact identity mismatch'
    print(a.command + ': exact 13-file release set verified')


if __name__ == '__main__':
    main()
