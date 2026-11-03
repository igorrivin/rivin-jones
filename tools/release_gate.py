#!/usr/bin/env python3
"""Check source publication requirements; hosted CI separately gates binaries."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    state = json.loads((ROOT / 'docs/release-state.json').read_text())
    pending = [item['description'] for item in state['source_requirements']
               if item['status'] != 'resolved']
    metadata = (ROOT / 'pyproject.toml').read_text()
    if 'Private :: Do Not Upload' in metadata:
        pending.append('Distribution metadata still marks the project private.')
    if 'license = "GPL-3.0-or-later"' not in metadata:
        pending.append('Project license metadata is missing.')
    for name in ('LICENSE', 'THIRD_PARTY_NOTICES.md', 'docs/REBUILDING.md',
                 'tools/prepare_vendor.py'):
        if not (ROOT / name).is_file():
            pending.append('Missing required source file: ' + name)
    if pending:
        print('\n'.join(pending), file=sys.stderr)
        return 1
    print('Source publication requirements passed. Binary publication also requires')
    print('the complete source archive and successful hosted tests of the exact artifacts.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
