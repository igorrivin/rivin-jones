"""Reusable contractions of the reduced colored R network, q=A^-4.

Only topology and contraction order are retained by a plan. Coefficients use
our division-free R matrices and highest-weight partial trace. Scalar execution
preserves the input field; optional NumPy batches explicitly use complex128.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product, islice

import cotengra

from .colored import _braid, _integer, _parameter, ColoredResourceLimitError
from .r_matrix import _Weights


@dataclass(frozen=True)
class _Crossing:
    sign: int
    legs: tuple
    indices: tuple
    trace_left: bool
    trace_right: bool


@dataclass(frozen=True)
class _Step:
    left: int
    right: int
    common_left: tuple
    common_right: tuple
    free_left: tuple
    free_right: tuple


def _topology(word, strands):
    pos = list(range(strands))
    raw = []
    for t, g in enumerate(word):
        p = abs(g)-1
        a, b = pos[p:p+2]
        c, d = strands+2*t, strands+2*t+1
        raw.append((1 if g > 0 else -1, (a, b, c, d),
                    p != 0 and a == p, b == p+1))
        pos[p:p+2] = c, d
    # Cut strand zero and fix both endpoints to highest weight (state zero).
    rename = dict(zip(pos, range(strands)))
    result = []
    for sign, legs, ta, tb in raw:
        legs = tuple(rename.get(x, x) for x in legs)
        legs = tuple(None if x == 0 else x for x in legs)
        counts = Counter(legs)
        # An index appearing twice in one leaf is traced inside that leaf.
        indices = tuple(x for x in legs if x is not None and counts[x] == 1)
        result.append(_Crossing(sign, legs, indices, ta, tb))
    counts = Counter(x for c in result for x in c.indices)
    if any(n != 2 for n in counts.values()):
        raise RuntimeError('colored network has an unclosed wire')
    return tuple(result)


def _schedule(inputs, path, color, sliced=()):
    """Compile binary joins; bound live dense entries including reshape copies."""
    sliced = frozenset(sliced)
    live = [tuple(x for x in inds if x not in sliced) for inds in inputs]
    sizes = [color**len(inds) for inds in live]
    total = sum(sizes)
    peak = total
    operations = 0
    steps = []
    if len(path) != max(0, len(live)-1):
        raise ValueError('contraction path must merge every tensor')
    for pair in path:
        if len(pair) != 2:
            raise ValueError('only binary contraction paths are supported')
        i, j = pair
        if not (0 <= i < j < len(live)):
            raise ValueError('invalid contraction path indices')
        a, b = live[i], live[j]
        common = tuple(x for x in a if x in b)
        ca, cb = tuple(a.index(x) for x in common), tuple(b.index(x) for x in common)
        fa = tuple(k for k, x in enumerate(a) if x not in common)
        fb = tuple(k for k, x in enumerate(b) if x not in common)
        output = tuple(a[k] for k in fa)+tuple(b[k] for k in fb)
        out_size = color**len(output)
        # Existing inputs remain live while a result and contiguous matrix
        # reshapes are allocated. Count both possible reshape copies.
        peak = max(peak, total+sizes[i]+sizes[j]+out_size)
        operations += out_size*color**len(common)
        total += out_size-sizes[i]-sizes[j]
        steps.append(_Step(i, j, ca, cb, fa, fb))
        for k in (j, i):
            live.pop(k)
            sizes.pop(k)
        live.append(output)
        sizes.append(out_size)
    if live and live[0]:
        raise RuntimeError('contraction leaves open indices')
    return tuple(steps), peak, operations


def _terms(crossing, color, fixed):
    """Generate only admissible conserved-charge R entries for one leaf."""
    n = color-1
    legs = crossing.legs
    limits = tuple(0 if x is None else fixed.get(x) for x in legs)
    irange = range(color) if limits[0] is None else (limits[0],)
    jrange = range(color) if limits[1] is None else (limits[1],)
    indices = tuple(x for x in crossing.indices if x not in fixed)
    positions = tuple(legs.index(x) for x in indices)
    for i in irange:
        for j in jrange:
            if crossing.sign > 0:
                low, high = 0, min(i, n-j)
                if limits[2] is not None:
                    low, high = max(low, limits[2]-j), min(high, limits[2]-j)
                if limits[3] is not None:
                    low, high = max(low, i-limits[3]), min(high, i-limits[3])
            else:
                low, high = 0, min(j, n-i)
                if limits[2] is not None:
                    low, high = max(low, j-limits[2]), min(high, j-limits[2])
                if limits[3] is not None:
                    low, high = max(low, limits[3]-i), min(high, limits[3]-i)
            for k in range(low, high+1):
                u, v = (j+k, i-k) if crossing.sign > 0 else (j-k, i+k)
                states = (i, j, u, v)
                if any(legs[a] is not None and legs[a] == legs[b] and states[a] != states[b]
                       for a in range(4) for b in range(a)):
                    continue
                exponent = 2*((n-2*i)*crossing.trace_left+(n-2*j)*crossing.trace_right)
                yield tuple(states[p] for p in positions), i, j, k, exponent


def _scalar_leaf(crossing, color, fixed, weights):
    result = {}
    for key, i, j, k, exponent in _terms(crossing, color, fixed):
        value = weights.coefficient(crossing.sign, i, j, k)*weights.power(exponent)
        if key in result:
            value = result[key]+value
        if value == weights.zero:
            result.pop(key, None)
        else:
            result[key] = value
    return result


def _join(left, right, step, zero):
    groups = defaultdict(list)
    for key, value in right.items():
        common = tuple(key[p] for p in step.common_right)
        groups[common].append((tuple(key[p] for p in step.free_right), value))
    result = {}
    for key, value in left.items():
        common = tuple(key[p] for p in step.common_left)
        prefix = tuple(key[p] for p in step.free_left)
        for suffix, other in groups.get(common, ()):
            output = prefix+suffix
            term = value*other
            if output in result:
                term = result[output]+term
            if term == zero:
                result.pop(output, None)
            else:
                result[output] = term
    return result


@dataclass(frozen=True)
class ColoredJonesPlan:
    """Immutable topology/path; construct with plan_colored_jones.

    max_entries bounds scalar slots in tensors and coefficient caches, using
    conservative dense bounds even for sparse execution. It does not bound
    Python overhead or the bit size of exact scalars. NumPy execution also
    applies a byte budget to estimated array payloads. No result is cached.
    """
    word: tuple
    strands: int
    color: int
    crossings: tuple
    path: tuple
    sliced_indices: tuple
    steps: tuple
    peak_entries: int
    estimated_operations: int
    cache_limit: int

    @property
    def slices(self):
        return self.color**len(self.sliced_indices)

    def info(self):
        return dict(algorithm='tensor', color=self.color, strands=self.strands,
                    crossings=len(self.word), slices=self.slices,
                    sliced_indices=self.sliced_indices, peak_entries=self.peak_entries,
                    estimated_operations=self.estimated_operations)

    def _framing(self, A, normalize):
        if not isinstance(normalize, bool):
            raise TypeError('normalize must be bool')
        n = self.color-1
        writhe = sum(1 if g > 0 else -1 for g in self.word)
        return A**(-writhe*n*(n+2)) if normalize else A**0*((-1)**(writhe*n % 2))

    def evaluate(self, A, *, normalize=True, return_info=False):
        """Evaluate in A's field (Fraction, exact field, Arb, or float/complex)."""
        A = _parameter(A)
        if not isinstance(return_info, bool):
            raise TypeError('return_info must be bool')
        framing = self._framing(A, normalize)
        if not self.crossings:
            result = A**0
            return (result, self.info()) if return_info else result
        weights = _Weights(A, self.color, self.cache_limit)
        total = weights.zero
        for values in product(range(self.color), repeat=len(self.sliced_indices)):
            fixed = dict(zip(self.sliced_indices, values))
            live = [_scalar_leaf(c, self.color, fixed, weights) for c in self.crossings]
            for step in self.steps:
                result = _join(live[step.left], live[step.right], step, weights.zero)
                live.pop(step.right)
                live.pop(step.left)
                live.append(result)
            total = total+live[0].get((), weights.zero)
        result = total*framing
        return (result, self.info()) if return_info else result

    def evaluate_many(self, As, *, backend='scalar', normalize=True, batch_size=128,
                      max_bytes=64*1024*1024, max_points=100000, return_info=False):
        """Reuse the plan for A values; NumPy batches are complex128 approximations.

        scalar returns a list preserving each input field. numpy returns a 1D
        complex128 array and chunks to fit max_bytes (array payload estimate,
        excluding interpreter/BLAS overhead). Exact/Arb inputs should use scalar.
        """
        batch_size = _integer(batch_size, 'batch_size', 1)
        max_bytes = _integer(max_bytes, 'max_bytes', 1)
        max_points = _integer(max_points, 'max_points', 1)
        if backend not in ('scalar', 'numpy'):
            raise ValueError('backend must be scalar or numpy')
        if not isinstance(normalize, bool) or not isinstance(return_info, bool):
            raise TypeError('normalize and return_info must be bool')
        values = list(islice(iter(As), max_points+1))
        if len(values) > max_points:
            raise ColoredResourceLimitError('evaluation batch exceeds max_points')
        info = dict(self.info(), backend=backend, points=len(values))
        if backend == 'scalar':
            result = [self.evaluate(A, normalize=normalize) for A in values]
        else:
            try:
                import numpy as np
            except ImportError as exc:
                raise ImportError('NumPy batches require the rivin-jones[tensor] extra') from exc
            # Inputs, outputs and per-chunk coefficient/contraction buffers.
            available = max_bytes-32*len(values)
            chunk = min(batch_size, available//(16*max(1, self.peak_entries)))
            if values and chunk < 1:
                raise ColoredResourceLimitError('max_bytes cannot hold one planned NumPy evaluation')
            if not values:
                chunk = 1
            points = np.asarray([complex(_parameter(A)) for A in values], dtype=np.complex128)
            if not np.all(np.isfinite(points)) or np.any(points == 0):
                raise ValueError('NumPy parameters must be finite and nonzero in complex128')
            result = np.empty(len(points), dtype=np.complex128)
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                for start in range(0, len(points), chunk):
                    result[start:start+chunk] = self._numpy(points[start:start+chunk], normalize, np)
            if not np.all(np.isfinite(result)):
                raise ArithmeticError('nonfinite tensor result; use exact or Arb scalar evaluation')
            info.update(batch_size=chunk, estimated_peak_bytes=32*len(values)+16*self.peak_entries*min(chunk, len(values)))
        return (result, info) if return_info else result

    def _numpy(self, points, normalize, np):
        count, N = len(points), self.color
        total = np.zeros(count, dtype=np.complex128)
        if not self.crossings:
            return total+1
        weights = _Weights(points, self.color, self.cache_limit)
        for values in product(range(N), repeat=len(self.sliced_indices)):
            fixed = dict(zip(self.sliced_indices, values))
            live = []
            for c in self.crossings:
                indices = tuple(x for x in c.indices if x not in fixed)
                array = np.zeros((count,)+(N,)*len(indices), dtype=np.complex128)
                for key, i, j, k, exponent in _terms(c, N, fixed):
                    array[(slice(None),)+key] += weights.coefficient(c.sign, i, j, k)*weights.power(exponent)
                live.append(array)
            del array
            for step in self.steps:
                left, right = live[step.left], live[step.right]
                # Batch axis is retained, never summed with a wire index.
                la = (0,)+tuple(p+1 for p in step.free_left+step.common_left)
                ra = (0,)+tuple(p+1 for p in step.common_right+step.free_right)
                a, b, k = N**len(step.free_left), N**len(step.free_right), N**len(step.common_left)
                out = (left.transpose(la).reshape(count, a, k) @
                       right.transpose(ra).reshape(count, k, b))
                out = out.reshape((count,)+(N,)*(len(step.free_left)+len(step.free_right)))
                live.pop(step.right)
                live.pop(step.left)
                live.append(out)
                del left, right
            total += live[0].reshape(count)
        return total*self._framing(points, normalize)


def plan_colored_jones(word, *, strands, color, optimize='greedy',
                       max_color=256, max_crossings=1000, max_entries=1000000,
                       max_slices=4096, max_operations=100000000, cache_limit=20000):
    """Compile an explicit experimental reduced colored contraction plan.

    optimize accepts Cotengra strategies/optimizers or a binary linear path;
    'sequential' supplies a no-search comparison. Slicing fixes selected wire
    indices and sums every assignment, without truncation or division. Planning
    rejects estimated work above max_operations before constructing coefficients.
    Hyperoptimizer search budgets are controlled by the supplied optimizer.
    """
    word, strands, components = _braid(word, strands)
    color = _integer(color, 'color', 1)
    max_color = _integer(max_color, 'max_color', 1)
    max_crossings = _integer(max_crossings, 'max_crossings')
    max_entries = _integer(max_entries, 'max_entries', 1)
    max_slices = _integer(max_slices, 'max_slices', 1)
    max_operations = _integer(max_operations, 'max_operations', 1)
    cache_limit = _integer(cache_limit, 'cache_limit')
    if components != 1:
        raise ValueError('colored Jones currently supports knots, not multi-component links')
    if color > max_color or len(word) > max_crossings:
        raise ColoredResourceLimitError('colored plan exceeds color or crossing limit')
    crossings = _topology(word, strands) if word and color > 1 else ()
    inputs = tuple(c.indices for c in crossings)
    indices = sorted({x for inds in inputs for x in inds})
    if len(inputs) <= 1:
        path = ()
    elif isinstance(optimize, str) and optimize == 'sequential':
        path = ((0, 1),)+tuple((0, n-1) for n in range(len(inputs)-1, 1, -1))
    else:
        tree = cotengra.array_contract_tree(inputs, output=(), size_dict={x: color for x in indices},
                                            optimize=optimize, canonicalize=False)
        if tree.sliced_inds:
            raise ValueError('pass an unsliced optimizer; the plan applies its own checked slicing')
        path = tuple(tuple(sorted(pair)) for pair in tree.get_path())
    # Gaussian binomials/factors plus bounded power/coefficient caches and
    # scratch vectors. Scalar object sizes and Python container overhead vary.
    reserve = (color*(color+1)+min(cache_limit, 2*color**3)
               +min(cache_limit, 12*(color-1)**2+1)+32)
    if not crossings:
        reserve = 1
    minimum_peak = len(crossings)+(3 if len(crossings) > 1 else 0)
    if reserve+minimum_peak > max_entries:
        raise ColoredResourceLimitError('max_entries cannot hold coefficient caches and scalar leaves; reduce cache_limit or raise max_entries')
    sliced = ()
    while True:
        steps, peak, cost = _schedule(inputs, path, color, sliced)
        if peak+reserve <= max_entries:
            break
        if color**(len(sliced)+1) > max_slices or len(sliced) == len(indices):
            raise ColoredResourceLimitError('tensor memory bound needs more slices; raise max_entries/max_slices or reduce cache_limit')
        # Minimize live memory first, arithmetic second. Fixed deterministic
        # tie-breaks make a plan reproducible for a supplied contraction path.
        choices = []
        for x in indices:
            if x not in sliced:
                _, p, c = _schedule(inputs, path, color, sliced+(x,))
                choices.append((p, c, x))
        sliced += (min(choices)[2],)
    # Enumerating local sparse supports takes at most N^3 trials per crossing.
    operations = (cost+len(crossings)*color**3)*color**len(sliced)
    if operations > max_operations:
        raise ColoredResourceLimitError(f'planned work {operations} exceeds max_operations={max_operations}')
    return ColoredJonesPlan(tuple(word), strands, color, crossings, path, sliced,
                            steps, peak+reserve, operations, cache_limit)
