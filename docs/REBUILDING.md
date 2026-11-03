# Rebuilding the native extension

Each release's source distribution, `rivin_jones-VERSION.tar.gz`, is the
corresponding source for its wheels. It is available alongside those wheels
on PyPI and in the [GitHub release](https://github.com/igorrivin/rivin-jones/releases).
It includes the Python and Rust application, every crate in Cargo.lock, the
GMP source bundled by gmp-mpfr-sys, license notices, tests and build scripts.
Ordinary Git checkouts omit `vendor/`; the release source archive includes it.

The build tools are Python 3.10–3.13, Rust/Cargo, Maturin 1.16.x, a C compiler,
GNU make, m4 and the usual POSIX shell/build utilities. On macOS, install the
Xcode Command Line Tools. On Linux, the standard C development tools are
needed. System compilers and the Python interpreter are not bundled. Release
CI records the actual compiler and platform versions with its artifacts.
Windows and Python 3.14+ are not supported by this initial release.

After extracting the source archive and entering its directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install 'maturin>=1.16,<2'
maturin build --release --locked --offline --interpreter python --out dist
python -m pip install dist/*.whl
python -m unittest discover -s tests -v
python -m pip check
```

Cargo's `--offline` uses the included `vendor/` sources through the included
`.cargo/config.toml`. Installing build tools and the separately distributed
Cotengra/Autoray runtime dependencies still needs network access or local
wheels. The Cargo registry cache is not needed for the native rebuild.
For a build independent of existing Cargo/GMP outputs, use fresh directories:

```sh
export CARGO_HOME="$(mktemp -d)"
export CARGO_TARGET_DIR="$(mktemp -d)"
export GMP_MPFR_SYS_CACHE="$(mktemp -d)"
maturin build --release --locked --offline --interpreter python --out dist
```

To modify an LGPL component, edit its source in `vendor/`, then make it a
local path dependency so Cargo does not enforce its original registry-file
checksums. For example, append to the application's Cargo.toml:

```toml
[patch.crates-io]
gmp-mpfr-sys = { path = "vendor/gmp-mpfr-sys-1.7.1" }
rug = { path = "vendor/rug-1.30.0" }
```

Use only the overrides for components you change. GMP's C sources are under
`vendor/gmp-mpfr-sys-1.7.1/gmp-6.3.0-c/`. After editing them, update the lockfile
for the path overrides and rebuild with a fresh target and GMP cache:

```sh
cargo update --offline
export CARGO_TARGET_DIR="$(mktemp -d)"
export GMP_MPFR_SYS_CACHE="$(mktemp -d)"
maturin build --release --locked --offline --interpreter python --out dist
python -m pip install --force-reinstall --no-deps dist/*.whl
python -m unittest discover -s tests -v
```

The replacement wheel contains the relinked extension. The project's GPL
terms permit this rebuilding and modification; the LGPL components retain
their own terms. No installation keys or locked device are involved.

To prepare a complete source archive from a Git checkout, first run
`python tools/prepare_vendor.py`, then `python tools/build_sdist.py --out dist` and
`python tools/check_distribution.py --public dist/*.tar.gz`. Do not publish
a wheel without making its matching complete source archive available.
