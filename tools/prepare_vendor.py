#!/usr/bin/env python3
"""Prepare complete native dependency sources for a release source archive."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
CONFIG = '[source.crates-io]\nreplace-with = "vendored-sources"\n\n[source.vendored-sources]\ndirectory = "vendor"\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true',
                        help='use only the existing Cargo cache')
    args = parser.parse_args()
    if (ROOT / 'vendor').exists() or (ROOT / '.cargo/config.toml').exists():
        parser.error('Use a fresh source checkout; vendor/config already exists.')
    command = ['cargo', 'vendor', '--locked', '--versioned-dirs']
    if args.offline:
        command.append('--offline')
    command.append('vendor')
    subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    (ROOT / '.cargo').mkdir(exist_ok=True)
    (ROOT / '.cargo/config.toml').write_text(CONFIG)
    print('Prepared locked sources and a portable Cargo configuration.')


if __name__ == '__main__':
    main()
