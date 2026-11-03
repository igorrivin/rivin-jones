#!/usr/bin/env python3
"""Build the complete source archive, retaining executable native build scripts."""
import argparse
import io
from pathlib import Path, PurePosixPath
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=ROOT/'dist')
    p.add_argument('--maturin', default='maturin')
    a = p.parse_args()
    out = a.out.resolve()
    assert (ROOT/'vendor').is_dir(), 'Run tools/prepare_vendor.py first'
    subprocess.run([a.maturin, 'sdist', '--out', str(out)], cwd=ROOT, check=True)
    archives = list(out.glob('rivin_jones-*.tar.gz'))
    assert len(archives) == 1, 'Use an output directory containing one release version'
    archive = archives[0]
    temporary = archive.with_suffix('.tmp')
    executable = 0
    # Maturin's extra-file glob currently writes regular files with mode 0644.
    # Native configure/install scripts need their original executable bits.
    with tarfile.open(archive, 'r:gz') as source, tarfile.open(temporary, 'w:gz') as target:
        for member in source:
            path = PurePosixPath(member.name)
            assert not path.is_absolute() and '..' not in path.parts
            assert not member.issym() and not member.islnk()
            original = ROOT.joinpath(*path.parts[1:])
            if member.isfile():
                if original.is_file() and original.stat().st_mode & 0o111:
                    member.mode = 0o755
                    executable += 1
                else:
                    member.mode = 0o644
                target.addfile(member, io.BytesIO(source.extractfile(member).read()))
            else:
                target.addfile(member)
    temporary.replace(archive)
    print('Preserved executable modes for', executable, 'source files in', archive)


if __name__ == '__main__':
    main()
