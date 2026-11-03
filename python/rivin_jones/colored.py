"""Reduced N-dimensional colored Jones of knot closures, q=t=A^-4.

The annular closure of f_n is S_n, where S_0=1, S_1=x and
S_(n+1)=x*S_n-S_(n-1). For n=N-1, evaluate S_n on blackboard cables,
divide EXACTLY by Delta_n=(-1)^n*sum_j A^(2*n-4*j), and correct the
framing by ((-1)^n*A^(n*(n+2)))^(-w). Specialize only afterwards.

The cable route uses native ordinary Jones calculations and exact division.
The direct route uses sparse colored R matrices and an open-strand quantum
trace, with native integer Laurent arithmetic or the supplied scalar field.
Neither route calls Sage for a polynomial or divides numerical zeros. The
explicit projector routine is separate and retains its coefficient field.
This module supports knots, not multi-component colored links.

References: Dasbach--Lin, Compositio Math. 142 (2006), Section 3;
https://doi.org/10.1112/S0010437X06002296 . The bracket here assigns delta to
one circle and 1 to the empty diagram, unlike their reduced bracket notation.
"""

from fractions import Fraction
from numbers import Integral
import math
import operator
import sys

from . import TLDiagram, jones_polynomial

DEFAULT_MAX_CABLED_STRANDS = 12
DEFAULT_MAX_CABLED_CROSSINGS = 500
_METHODS = ("auto", "frontier", "frontier-gmp", "greedy", "sequential")


class ColoredResourceLimitError(ValueError):
    """A cable, projector, polynomial or direct contraction exceeds its guard."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or getattr(getattr(value, 'dtype', None), 'kind', None) == 'b':
        raise TypeError(f"{name} must be an integer, not bool")
    value = operator.index(value)
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _braid(word, strands):
    strands = _integer(strands, "strands", 1)
    if strands > 128:
        raise ValueError("strands must be at most 128")
    generators = []
    permutation = list(range(strands))
    for g in word:
        if isinstance(g, bool) or getattr(getattr(g, 'dtype', None), 'kind', None) == 'b':
            raise TypeError("braid generators must be integers, not bool")
        g = operator.index(g)
        if not 1 <= abs(g) < strands:
            raise ValueError("braid generators must satisfy 1 <= abs(g) < strands")
        generators.append(g)
        i = abs(g)-1
        permutation[i], permutation[i+1] = permutation[i+1], permutation[i]
    unseen = set(range(strands))
    components = 0
    while unseen:
        components += 1
        i = next(iter(unseen))
        while i in unseen:
            unseen.remove(i)
            i = permutation[i]
    return generators, strands, components


def _cable(word, strands, width):
    result = []
    for g in word:
        base = (abs(g)-1)*width
        block = [base+k+j+1 for k in range(width-1, -1, -1) for j in range(width)]
        result.extend(block if g > 0 else [-h for h in reversed(block)])
    return result, strands*width


def cable_braid(braid_word, n_strands, cable_width, *, max_strands=128,
                max_crossings=100000):
    """Blackboard cable; a negative bundle crossing is the positive block inverse.

    Every original crossing becomes width**2 crossings. This swaps two bundles;
    it is not a full twist of the two bundles.
    """
    word, strands, _ = _braid(braid_word, n_strands)
    width = _integer(cable_width, "cable_width", 1)
    limit = _integer(max_strands, "max_strands", 1)
    crossings = _integer(max_crossings, "max_crossings")
    if strands*width > min(limit, 128) or len(word)*width*width > crossings:
        raise ColoredResourceLimitError("cabled braid exceeds the declared size limits")
    return _cable(word, strands, width)


def _add_scaled(target, poly, scale=1, shift=0):
    for exponent, coefficient in poly.items():
        key = exponent+shift
        value = target.get(key, 0)+scale*coefficient
        if value:
            target[key] = value
        else:
            target.pop(key, None)


def _shift(poly, exponent, coefficient=1):
    return {e+exponent: coefficient*c for e, c in poly.items() if c}


def _divide_exact(numerator, denominator):
    """Long division in Z[A,A^-1]; never discard a nonzero remainder."""
    if not denominator:
        raise ZeroDivisionError("zero Laurent polynomial divisor")
    if not numerator:
        return {}
    low, dlow = min(numerator), min(denominator)
    remainder = _shift(numerator, -low)
    divisor = _shift(denominator, -dlow)
    degree = max(divisor)
    answer = {}
    while remainder and max(remainder) >= degree:
        high = max(remainder)
        coefficient, residual = divmod(remainder[high], divisor[degree])
        if residual:
            raise ArithmeticError("colored closure division is not integral")
        exponent = high-degree
        answer[exponent+low-dlow] = coefficient
        _add_scaled(remainder, divisor, -coefficient, exponent)
    if remainder:
        raise ArithmeticError("colored closure has a nonzero quantum-dimension remainder")
    return answer


def _cabling_polynomial(word, *, strands, color, normalize=True,
                             method="auto", max_cabled_strands=DEFAULT_MAX_CABLED_STRANDS,
                             max_cabled_crossings=DEFAULT_MAX_CABLED_CROSSINGS):
    """Return exact integer coefficients {A_exponent: coefficient} for a knot.

    Color is the representation dimension N>=1; J_1=1 and J_2 is ordinary
    reduced Jones. normalize=False keeps blackboard framing, while still
    dividing by the single knot's quantum dimension. The default is the
    ambient-isotopy invariant. Limits are checked before cable allocation.

    Cabling grows rapidly with color; these bounds do not promise a runtime
    or memory bound. Increase them explicitly only for deliberate experiments.
    """
    color = _integer(color, "color", 1)
    limit = _integer(max_cabled_strands, "max_cabled_strands", 1)
    crossing_limit = _integer(max_cabled_crossings, "max_cabled_crossings")
    if limit > 128:
        raise ValueError("max_cabled_strands must be at most 128")
    if not isinstance(normalize, bool):
        raise TypeError("normalize must be bool")
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}")
    word, strands, components = _braid(word, strands)
    if components != 1:
        raise ValueError("colored Jones currently supports knots, not multi-component links")
    if color == 1 or not word:
        return {0: 1}
    n = color-1
    if strands*n > limit or len(word)*n*n > crossing_limit:
        raise ColoredResourceLimitError(
            f"cable requires {strands*n} strands and {len(word)*n*n} crossings; "
            f"limits are {limit} and {crossing_limit}")
    if method in ("frontier", "frontier-gmp") and strands*n > 12:
        raise ColoredResourceLimitError("the explicit frontier method supports at most 12 cabled strands")
    writhe = sum(1 if g > 0 else -1 for g in word)
    if color == 2:
        ordinary = jones_polynomial(word, strands=strands, method=method)
        return ordinary if normalize else _shift(ordinary, 3*writhe, (-1)**(writhe % 2))

    closure = {}
    for j in range(n//2+1):
        width = n-2*j
        weight = (-1)**j*math.comb(n-j, j)
        if width == 0:
            # Empty cable has unreduced bracket 1, not the one-circle value.
            _add_scaled(closure, {0: 1}, weight)
            continue
        cable, cable_strands = _cable(word, strands, width)
        ordinary = jones_polynomial(cable, strands=cable_strands, method=method)
        cable_writhe = writhe*width*width
        # Undo ordinary reduction and framing: B(D^width)=delta*(-A^3)^w_c*J_2.
        scale = -weight*((-1)**(cable_writhe % 2))
        _add_scaled(closure, ordinary, scale, 3*cable_writhe+2)
        _add_scaled(closure, ordinary, scale, 3*cable_writhe-2)
    dimension = {2*n-4*j: (-1)**n for j in range(n+1)}
    reduced = _divide_exact(closure, dimension)
    if normalize:
        reduced = _shift(reduced, -writhe*n*(n+2), (-1)**((n*writhe) % 2))
        if any(e % 4 for e in reduced):
            raise ArithmeticError("normalized knot polynomial left the q=A^-4 exponent lattice")
    return reduced


def _direct_options(word, strands, color, normalize, max_color, max_rows,
                    max_states, max_transitions, cache_limit):
    color = _integer(color, 'color', 1)
    max_color = _integer(max_color, 'max_color', 1)
    if not isinstance(normalize, bool):
        raise TypeError('normalize must be bool')
    limits = dict(max_rows=_integer(max_rows, 'max_rows', 1),
                  max_states=_integer(max_states, 'max_states', 1),
                  max_transitions=_integer(max_transitions, 'max_transitions', 1),
                  cache_limit=_integer(cache_limit, 'cache_limit'))
    word, strands, components = _braid(word, strands)
    if components != 1:
        raise ValueError('colored Jones currently supports knots, not multi-component links')
    if word and color > max_color:
        raise ColoredResourceLimitError(f'color {color} exceeds direct limit {max_color}')
    return word, strands, color, limits


def colored_jones_polynomial(word, *, strands, color, normalize=True, algorithm='auto',
                             method='auto', max_cabled_strands=DEFAULT_MAX_CABLED_STRANDS,
                             max_cabled_crossings=DEFAULT_MAX_CABLED_CROSSINGS,
                             max_color=32, max_rows=1000000, max_states=10000,
                             max_transitions=2000000, cache_limit=20000,
                             max_polynomial_span=200000):
    """Exact reduced colored knot polynomial using cables or sparse R matrices.

    auto keeps ordinary color 2 and small cables (<=4 strands and <=80
    crossings), and uses the direct representation otherwise. An explicit ordinary contraction
    method selects cabling unless algorithm='r-matrix' was explicitly chosen.
    Cable limits apply to cables; direct limits apply to R contraction.
    """
    if algorithm not in ('auto', 'cabling', 'r-matrix'):
        raise ValueError('algorithm must be auto, cabling or r-matrix')
    if method not in _METHODS:
        raise ValueError(f'method must be one of {_METHODS}')
    cable_limit = _integer(max_cabled_strands, 'max_cabled_strands', 1)
    if cable_limit > 128:
        raise ValueError('max_cabled_strands must be at most 128')
    crossing_limit = _integer(max_cabled_crossings, 'max_cabled_crossings')
    max_color = _integer(max_color, 'max_color', 1)
    span_limit = _integer(max_polynomial_span, 'max_polynomial_span', 1)
    # Validate input without applying the direct color bound to the cable path.
    word, strands, color, limits = _direct_options(word, strands, color, normalize,
        max(_integer(color, 'color', 1), max_color),
        max_rows, max_states, max_transitions, cache_limit)
    n = color-1
    if color == 1 or not word:
        return {0: 1}
    if algorithm == 'auto':
        algorithm = 'cabling' if color == 2 or method != 'auto' or (strands*n <= 4 and len(word)*n*n <= 80) else 'r-matrix'
    if algorithm == 'cabling':
        return _cabling_polynomial(word, strands=strands, color=color, normalize=normalize,
            method=method, max_cabled_strands=cable_limit, max_cabled_crossings=crossing_limit)
    if color > max_color:
        raise ColoredResourceLimitError(f'color {color} exceeds direct limit {max_color}')
    bound = 12*len(word)*n*n+4*(strands-1)*n
    if bound > span_limit:
        raise ColoredResourceLimitError(f'polynomial span bound {bound} exceeds {span_limit}')
    from .r_matrix import _contract
    from . import _tl_tensor
    value, _ = _contract(word, strands, color, _tl_tensor.PolyWrapper(1, [1]),
                          normalize=normalize, **limits)
    result = {e: c for c, e in value.terms}
    if normalize and any(e % 4 for e in result):
        raise ArithmeticError('normalized knot polynomial left the q=A^-4 exponent lattice')
    return result


def _parameter(A):
    if isinstance(A, bool) or getattr(getattr(A, 'dtype', None), 'kind', None) == 'b':
        raise TypeError("A must be a nonzero scalar, not bool")
    if isinstance(A, Integral):
        A = Fraction(int(A))
    if A == 0:
        raise ValueError("A must be nonzero")
    if isinstance(A, (float, complex)) or getattr(getattr(A, 'dtype', None), 'kind', None) in ('f', 'c'):
        z = complex(A)
        if not math.isfinite(z.real) or not math.isfinite(z.imag):
            raise ValueError("A must be finite")
    if hasattr(A, "contains_zero") and A.contains_zero():
        raise ValueError("the A enclosure must exclude zero")
    return A


def evaluate_polynomial(polynomial, A):
    """Evaluate a Laurent coefficient map, preserving A's scalar arithmetic.

    Fractions/exact fields give exact values; Arb fields give enclosures;
    Python float/complex gives an approximation without an error certificate.
    """
    A = _parameter(A)
    return sum((c*A**e for e, c in polynomial.items()), A**0*0)


def colored_jones_evaluate(word, *, strands, color, A, algorithm='auto',
                           return_info=False, normalize=True, max_color=256,
                           max_rows=1000000, max_states=10000, max_transitions=2000000,
                           cache_limit=20000, **options):
    """Direct reduced open-strand trace, regular at roots, in A's scalar field.

    algorithm='polynomial' computes an exact polynomial first; 'cabling' forces
    the earlier cable reference. auto and r-matrix evaluate directly without
    expanding either cables or a complete polynomial. return_info=True returns
    (value, work_counters). No cross-call result or contraction cache is used.
    tensor compiles a bounded contraction plan; use plan_colored_jones directly
    to reuse it. Tensor options are optimize, max_entries, max_slices,
    max_operations and max_crossings; row/state/transition limits apply only
    to the sequential R-matrix route.
    """
    A = _parameter(A)
    if not isinstance(return_info, bool):
        raise TypeError('return_info must be bool')
    if algorithm in ('polynomial', 'cabling'):
        polynomial = colored_jones_polynomial(word, strands=strands, color=color,
            algorithm='auto' if algorithm == 'polynomial' else 'cabling', normalize=normalize,
            max_color=max_color, max_rows=max_rows, max_states=max_states,
            max_transitions=max_transitions, cache_limit=cache_limit, **options)
        value = evaluate_polynomial(polynomial, A)
        return (value, dict(algorithm=algorithm)) if return_info else value
    if algorithm == 'tensor':
        from .tensor import plan_colored_jones
        plan = plan_colored_jones(word, strands=strands, color=color,
                                 max_color=max_color, cache_limit=cache_limit, **options)
        return plan.evaluate(A, normalize=normalize, return_info=return_info)
    if algorithm not in ('auto', 'r-matrix'):
        raise ValueError('evaluation algorithm must be auto, r-matrix, tensor, polynomial or cabling')
    if options:
        raise TypeError('direct evaluation does not accept cable or polynomial options: '+', '.join(options))
    word, strands, color, limits = _direct_options(word, strands, color, normalize, max_color,
        max_rows, max_states, max_transitions, cache_limit)
    if color == 1 or not word:
        value, info = A**0, dict(algorithm='trivial', color=color, strands=strands)
    else:
        from .r_matrix import _contract
        value, info = _contract(word, strands, color, A, normalize=normalize, **limits)
    return (value, info) if return_info else value


def colored_jones_cabling(braid_word, n_strands, N, A, normalize=True, **options):
    """Compatibility entry point; now computes the actual reduced colored invariant."""
    return colored_jones_evaluate(braid_word, strands=n_strands, color=N, A=A,
                                  normalize=normalize, algorithm='cabling', **options)


def colored_jones_numeric(braid_word, n_strands, N, A, normalize=True, **options):
    """Compatibility entry point for the regular direct colored evaluation."""
    return colored_jones_evaluate(braid_word, strands=n_strands, color=N, A=A,
                                  normalize=normalize, **options)


def quantum_integer(n, A):
    """[n] for x=-A^2: a finite sum, regular also at x=+/-1; [2]=delta."""
    n = _integer(n, "n")
    A = _parameter(A)
    x = -A**2
    return sum((x**(n-1-2*j) for j in range(n)), A**0*0)


def _tl_product(left, right, delta):
    result = {}
    for idx, coefficient in left.items():
        diagram = TLDiagram(list(idx))
        for other_idx, other_coefficient in right.items():
            product = diagram @ TLDiagram(list(other_idx))
            key = tuple(product.idx)
            value = result.get(key, 0)+coefficient*other_coefficient*delta**product.loops
            if value == 0:
                result.pop(key, None)
            else:
                result[key] = value
    return result


def jones_wenzl_projector(n, A, *, max_strands=6):
    """f_n as {pairing_tuple: coefficient}, preserving exact/ball arithmetic.

    Recursion requires [2],...,[n] to be invertible. At the Volume Conjecture
    root for color N, f_(N-1) exists but its trace [N] vanishes. No coefficient
    is dropped using a magnitude threshold. This routine does not divide by
    the trace. A ball denominator containing zero or a numerically unresolved
    floating denominator is rejected.
    """
    n = _integer(n, "n")
    limit = _integer(max_strands, "max_strands")
    if n > min(limit, 128):
        raise ColoredResourceLimitError("projector exceeds the declared strand limit")
    A = _parameter(A)
    one = A**0
    projector = {(): one}
    delta = -A**2-A**(-2)
    for k in range(1, n+1):
        extended = {tuple((TLDiagram(list(idx)) | TLDiagram.eye(1)).idx): coefficient
                    for idx, coefficient in projector.items()}
        if k == 1:
            projector = extended
            continue
        denominator = quantum_integer(k, A)
        if denominator == 0 or (hasattr(denominator, "contains_zero") and denominator.contains_zero()):
            raise ZeroDivisionError(f"projector f_{n} needs nonzero [{k}]")
        if isinstance(denominator, (float, complex)):
            x = -A**2
            scale = sum(abs(x**(k-1-2*j)) for j in range(k))
            if not math.isfinite(abs(denominator)) or abs(denominator) <= 32*sys.float_info.epsilon*scale:
                raise ArithmeticError(f"floating arithmetic cannot resolve projector denominator [{k}]")
        ratio = quantum_integer(k-1, A)/denominator
        generator = {tuple(TLDiagram.gen(k, k-2).idx): one}
        product = _tl_product(_tl_product(extended, generator, delta), extended, delta)
        projector = dict(extended)
        for idx, coefficient in product.items():
            value = projector.get(idx, 0)-ratio*coefficient
            if value == 0:
                projector.pop(idx, None)
            else:
                projector[idx] = value
    return projector
