The first public release of rivin-jones computes exact ordinary and colored
Jones polynomials, evaluates colored knot invariants directly at roots of
unity, and supports reusable tensor contraction plans with exact slicing and
optional NumPy batches.

The recorded 35-case complete-polynomial benchmark was 5.61–3,364 times faster
than Sage 10.9 on the same macOS arm64 machine, with complete coefficient
agreement. These are finite-corpus measurements of the native candidate,
not a universal speed guarantee. New tensor-plan measurements compare against
our previous engine and are reported separately in docs/PLANNED_CONTRACTION.md.

Python 3.10–3.13 wheels are built and tested for Linux x86_64, macOS arm64 and
macOS x86_64. Sage is optional. The source archive includes the locked native
dependency sources, GMP and rebuild instructions. See docs/REBUILDING.md.

License: GPL-3.0-or-later, with retained third-party notices.

```sh
python -m pip install rivin-jones
# Optional complex128 batches:
python -m pip install 'rivin-jones[tensor]'
```

Polynomial outputs use integer powers of A, with q = A^-4. Color N denotes
representation dimension; the unknot is normalized to 1. Colored input must
close to a knot. Arb values enclose finite-color evaluations, not the error
to the conjectured limiting volume.
