# Third-party terms and dependencies

The combined rivin-jones distribution is GPL-3.0-or-later; see LICENSE.
Portions received under CC0 1.0 omit voluntary attribution at the original
maintainer's request. This does not remove any separate dependency's notices.

Cargo.lock pins the Rust dependencies. `licenses/dependencies.json` records
their license expressions and source hashes; `licenses/` retains their
license and notice texts. The complete release source archive also contains
the locked crates in `vendor/`, including their original notices and build
scripts. In particular:

- series 0.13.0 is GPL-3.0-or-later.
- Rug 1.30.0 and gmp-mpfr-sys 1.7.1 are LGPL-3.0-or-later. The latter builds
  bundled GMP 6.3.0, whose COPYING files and source notices are preserved.
- Other Rust dependencies use the MIT, Apache, Unicode and related terms
  listed in the inventory. These terms continue to apply to those components.

The native extension statically incorporates dependencies. For every release,
the corresponding source is the same-version `.tar.gz` source distribution
available beside the wheels on PyPI and at:

https://github.com/igorrivin/rivin-jones/releases

That archive contains the application and native dependency sources, Cargo
lockfile, compiler configuration and build/relink instructions. It can be
used to rebuild the extension with changed LGPL components; see
[docs/REBUILDING.md](docs/REBUILDING.md). No restriction on modification or
reverse engineering for debugging such modifications is added by this project.

Cotengra and Autoray are separately installed Python dependencies; they are
not copied into this wheel. NumPy and plotting packages are optional extras.
Sage is an optional independent oracle and provider of exact/Arb fields,
not a runtime dependency. No Sage Jones routine is called by the native engine.

Mathematical references and conventions appear in docs/MATHEMATICS.md.
