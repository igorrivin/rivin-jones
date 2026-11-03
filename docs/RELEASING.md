# Repository and distribution release

The repository is `igorrivin/rivin-jones`, distribution `rivin-jones`, import
`rivin_jones`. The project is GPL-3.0-or-later. Inherited CC0 portions omit
voluntary attribution at the original maintainer's request. Preserve the
separate dependency notices and complete corresponding source for each wheel.

## Source and local validation

`python tools/release_gate.py` checks source publication prerequisites.
The source repository can be public before hosted binary validation. A
source checkout excludes registry caches and build outputs. Do not add
private correspondence or historical workspace exports to this repository.

Prepare native dependency sources in a fresh checkout and build:

```sh
python tools/prepare_vendor.py
python tools/build_sdist.py --out dist
mkdir extracted
tar -xzf dist/*.tar.gz -C extracted --strip-components=1
(cd extracted && maturin build --release --locked --offline --out ../dist)
python tools/check_distribution.py --public dist/*
python -m twine check dist/*
python -m pip install dist/*.whl
python -m unittest discover -s tests -v
python -m pip check
```

The wheel is built from the generated source archive. `build_sdist.py` also
preserves executable modes on vendored native build scripts. The artifact
check verifies all vendored files against the Cargo checksums and requires
the GMP source and rebuild instructions. The source distribution is also
the complete corresponding-source bundle. See [REBUILDING.md](REBUILDING.md)
for an independent rebuild and a modified-library relink.

Set `RIVIN_JONES_TEST_SAGE=1` under Sage Python for additional exact and Arb
oracle tests. Sage is not a release runtime dependency.

## Hosted distributions

`distributions.yml` builds one complete source archive, then builds every
wheel from that archive, offline with respect to Cargo. Linux uses a
manylinux2014 container. The matrix covers Python 3.10, 3.11, 3.12 and 3.13
on Linux x86_64 and both macOS architectures. Each wheel is installed and
tested, including optional NumPy batches, before artifact collection.
The final job records the exact commit, run ID and 13 artifact hashes in
`SHA256SUMS.json` and uploads the immutable `release-distributions` artifact.

## PyPI and GitHub releases

Configure a [PyPI trusted publisher](https://docs.pypi.org/trusted-publishers/)
with owner `igorrivin`, repository `rivin-jones`, workflow `publish.yml`,
environment `pypi`. The first release uses an account-level pending publisher.
No long-lived PyPI token belongs in this repository.

Dispatch `publish.yml` on main with the successful distributions run ID.
It requires that the run succeeded at the current main commit, verifies
every artifact hash and source archive, then publishes those unchanged files.
It does not rebuild after validation. Only the isolated publishing job can
request a PyPI OIDC token.

After upload, the workflow downloads every public PyPI file and checks its
SHA-256, installs from the public index into a clean Linux/Python 3.12
environment, runs the suite and `pip check`, and creates the GitHub release
with the same files and hashes. The public verification report is retained
as an Actions artifact. Source publication, CI success, publisher setup and
public-index verification are separate results; report each accurately.
