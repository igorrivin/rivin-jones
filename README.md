# rivin-jones

Exact ordinary and colored Jones polynomials, with direct evaluation of
colored knot invariants at roots of unity.

Licensed under [GPL-3.0-or-later](LICENSE), with retained
[third-party notices](THIRD_PARTY_NOTICES.md). Release source archives include
the native dependency sources and [rebuild instructions](docs/REBUILDING.md).

The ordinary engine uses exact Temperley–Lieb frontier/tensor contraction.
The colored engine uses sparse R matrices and a reduced open-strand quantum
trace, with exact cabling retained as an independent calculation. No Sage
Jones routine is called during native computation.
An explicit planned-contraction route adds reusable paths, exact slicing,
and optional NumPy batches to the same colored R-matrix formulas.

## Installation

Python 3.10–3.13 is supported. Release wheels target Linux x86_64, macOS arm64
and macOS x86_64. Building from source needs Rust, a C compiler and the native
build tools described in [REBUILDING.md](docs/REBUILDING.md).

```sh
python -m pip install rivin-jones
# Optional NumPy batches:
python -m pip install 'rivin-jones[tensor]'
```

Cotengra is a runtime dependency. Sage and plotting packages are optional.
To build a Git checkout locally, use `python -m pip install .`.

## Exact polynomials

```python
from rivin_jones import jones_polynomial, colored_jones_polynomial

trefoil = jones_polynomial([1, 1, 1], strands=2)
assert trefoil == {-4: 1, -12: 1, -16: -1}

figure_eight = [1, -2, 1, -2]
colored = colored_jones_polynomial(figure_eight, strands=3, color=10)
```

Results are dictionaries `{A_exponent: integer_coefficient}`. The usual
Jones variable is `q=t=A^-4`. Color N is representation dimension, with
`J_N(unknot)=1` and `J_2` equal to the ordinary reduced Jones polynomial.
Braid width is explicit, including idle strands for ordinary links. Colored
input must close to one knot; multi-component colored links are rejected.

`algorithm='auto'` selects small cables or sparse R contraction. Force either
with `'cabling'` or `'r-matrix'`. Ordinary `method` options are `'auto'`,
`'frontier'`, `'frontier-gmp'`, `'greedy'` and `'sequential'`.

## Exact or interval evaluation

```python
from fractions import Fraction
from rivin_jones import colored_jones_evaluate

value = colored_jones_evaluate([1, 1, 1], strands=2, color=5,
                               A=Fraction(3, 2))
```

The scalar's arithmetic is preserved: rational/exact field inputs give exact
results; Arb inputs give enclosures; Python float/complex inputs give
approximations. Evaluation defaults to direct R contraction without the full
polynomial. `algorithm='polynomial'` constructs the polynomial first and
`'cabling'` selects the cable reference. `normalize=False` retains blackboard
framing with reduced normalization. `return_info=True` also returns work counts.

## Reusable colored tensor plans

Wider braids can benefit from a different contraction order. Build a plan
once and reuse it with new scalar values:

```python
from fractions import Fraction
from rivin_jones import plan_colored_jones

word = [1, -2, 3, -4] * 3
plan = plan_colored_jones(word, strands=5, color=4)
value = plan.evaluate(Fraction(3, 2))       # exact rational arithmetic
values = plan.evaluate_many([1j, -1j])    # preserves scalar arithmetic
print(plan.info())
```

`colored_jones_evaluate(..., algorithm='tensor')` builds and executes a plan
for a single call. The existing automatic evaluator remains the sparse braid
sweep, which wins on favorable narrow/high-color inputs. A plan holds only
topology and a path; it does not cache answers. The input word, width and color
are fixed by the plan; normalization can change at evaluation time.

Install the optional `tensor` extra to evaluate batches
with NumPy. This mode explicitly uses complex128 approximations:

```python
import numpy as np

As = np.exp(-0.2j + 1j*np.linspace(-0.02, 0.02, 32))
values, info = plan.evaluate_many(As, backend='numpy', return_info=True,
                                 batch_size=128, max_bytes=64*1024*1024)
```

Output arrays are complex128 even when the input array contains integers or
real values. Exact fields and Arb enclosures should use the scalar backend.
All methods use the open-strand reduction and remain defined at nonzero roots
of unity. Approximate arithmetic still requires conditioning/precision checks.

The planner accepts Cotengra strategies or an unsliced optimizer through
`optimize`; `'sequential'` is also available for comparison. `max_entries`
limits estimated live tensor and coefficient-cache slots. If necessary, the
plan fixes selected wire indices, evaluates each smaller network, and adds
**all** slice values. `max_slices` and `max_operations` bound that tradeoff.
NumPy batches additionally shrink their chunk size to fit the estimated array
payload budget `max_bytes`. These bounds exclude Python/BLAS overhead and the
variable bit sizes of exact scalar objects; they are not process RSS limits.
See [planned contraction details](docs/PLANNED_CONTRACTION.md).

For the Volume Conjecture, increase color N and evaluate at
`A_N=exp(-pi*i/(2*N))`. With Sage installed for Arb fields:

```python
from sage.all import RealBallField, ComplexBallField
from rivin_jones import colored_jones_evaluate

N = 100
R, C = RealBallField(512), ComplexBallField(512)
A = C(0, -R.pi()/(2*N)).exp()
value = colored_jones_evaluate([1, -2, 1, -2], strands=3, color=N, A=A)
assert abs(value) > 0
growth = 2*R.pi()/N*abs(value).log()
```

This encloses the finite-N value, not the error to the conjectured limiting
volume. The [CLI example](examples/volume_conjecture.py) records enclosures and
work counts. The [mathematical specification](docs/MATHEMATICS.md) states the
crossing matrices, inverse, trace, framing and pruning argument.

## Limits

Direct defaults are color <=32 for polynomials and <=256 for evaluations,
1,000,000 trace rows, 10,000 live states per row and 2,000,000 transitions.
The polynomial support-span bound is 200,000. Cables default to at most 12
strands and 500 crossings. Exceeding a guard raises
`rivin_jones.colored.ColoredResourceLimitError` without a partial answer.
Keyword arguments allow deliberate changes to the limits.

These are workload guards, not runtime or byte-memory guarantees. Direct
trace rows grow as N^(strands-1); the figure-eight presentation is unusually
favorable for pruning. Evidence at color 200 for that knot is not a promise
for arbitrary knots at that color.

## Evidence and development

The frozen candidate passed 405 direct-engine audit groups and 330 prior
colored regression groups. Its 35 repeated polynomial cases were faster
than Sage with exact coefficient equality; five completed root comparisons
were faster, and one additional Sage comparison remained incomplete.
See [performance and provenance](docs/PERFORMANCE.md). Those measurements
precede packaging; the renamed package is checked by coefficient replay and
fresh invariant tests, not advertised as a new timing sweep.

```sh
python -m unittest discover -s tests -v
python tools/check_distribution.py dist/*
python tools/release_gate.py
```

Portable tests require only the installed package. Set `RIVIN_JONES_TEST_SAGE=1`
under a Sage Python interpreter to enable the additional exact/Arb oracle
tests. The source gate checks licensing and notices; hosted CI gates binary publication.
[Release instructions](docs/RELEASING.md) separate local build, hosted CI,
repository publication and PyPI publication.
