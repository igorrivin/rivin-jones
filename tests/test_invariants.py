"""Portable exact checks; no Sage installation is needed."""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random
import subprocess
import sys
import unittest

import rivin_jones as jones
from rivin_jones import colored, r_matrix
import pairing_oracle as graph

DATA = Path(__file__).parent/'data'


def multiply(a, b):
    result = {}
    for e, c in a.items():
        for f, d in b.items(): result[e+f] = result.get(e+f, 0)+c*d
    return {e: c for e, c in result.items() if c}


def action(weights, word, initial):
    states = {initial: weights.one}
    for generator in word:
        p, sign = abs(generator)-1, 1 if generator > 0 else -1
        result = {}
        for state, amplitude in states.items():
            i, k = state[p:p+2]
            bound = min(i, weights.n-k) if sign > 0 else min(k, weights.n-i)
            for j in range(bound+1):
                pair = (k+j, i-j) if sign > 0 else (k-j, i+j)
                output = state[:p]+pair+state[p+2:]
                result[output] = result.get(output, weights.zero)+amplitude*weights.coefficient(sign, i, k, j)
        states = {s: a for s, a in result.items() if a != weights.zero}
    return states


def diagram_product(a, b, delta):
    result = {}
    for left, c in a.items():
        for right, d in b.items():
            pairing, loops = graph.compose(left, right)
            result[pairing] = result.get(pairing, 0)+c*d*delta**loops
    return {p: c for p, c in result.items() if c}


class Fixtures(unittest.TestCase):
    pass


def fixture_test(case, ordinary):
    def check(self):
        expected = {int(e): c for e, c in case['coefficients'].items()}
        if ordinary:
            actual = jones.jones_polynomial(case['word'], strands=case['strands'])
        else:
            options = dict(strands=case['strands'], color=case['color'])
            actual = jones.colored_jones_polynomial(case['word'], **options)
            self.assertEqual(jones.colored_jones_polynomial(case['word'], algorithm='r-matrix', **options), expected)
        self.assertEqual(actual, expected)
    return check


for kind in ('ordinary', 'colored'):
    for index, case in enumerate(json.loads((DATA/(kind+'.json')).read_text())):
        setattr(Fixtures, f'test_{kind}_{index:03d}_{case["name"]}', fixture_test(case, kind == 'ordinary'))


class Algebra(unittest.TestCase):
    def test_crossing_inverse_and_yang_baxter(self):
        A = Fraction(3, 2)
        for N in (2, 3, 4):
            w = r_matrix._Weights(A, N, 20000)
            for state in product(range(N), repeat=2):
                for word in ([1, -1], [-1, 1]):
                    self.assertEqual(action(w, word, state), {state: Fraction(1)})
            for state in product(range(N), repeat=3):
                for sign in (-1, 1):
                    self.assertEqual(action(w, [sign, 2*sign, sign], state), action(w, [2*sign, sign, 2*sign], state))

    def test_markov_mirror_framing_and_pruning(self):
        A = Fraction(4, 3)
        for word, strands in [([1]*3, 2), ([1, -2]*2, 3), ([1, -2, 1, 1, 2, 2], 3)]:
            for N in (2, 3, 4):
                def compute(w, s=strands, **kw):
                    return jones.colored_jones_polynomial(w, strands=s, color=N, algorithm='r-matrix', **kw)
                expected = compute(word)
                for changed, width in [(word+[strands], strands+1), (word+[-strands], strands+1),
                                       ([1]+word+[-1], strands), ([1, -1]+word, strands)]:
                    self.assertEqual(compute(changed, width), expected)
                self.assertEqual(compute([-g for g in word]), {-e: c for e, c in expected.items()})
                n, writhe = N-1, sum(1 if g > 0 else -1 for g in word)
                self.assertEqual(compute(word, normalize=False), {e+writhe*n*(n+2): (-1)**(writhe*n % 2)*c for e, c in expected.items()})
                value = colored.evaluate_polynomial(expected, A)
                for prune, cache_limit in [(True, 20000), (False, 20000), (True, 0)]:
                    actual, _ = r_matrix._contract(word, strands, N, A, prune=prune, cache_limit=cache_limit)
                    self.assertEqual(actual, value)

    def test_connected_sum(self):
        for N in (2, 3, 4, 5):
            a = jones.colored_jones_polynomial([1]*3, strands=2, color=N)
            b = {-e: c for e, c in a.items()}
            self.assertEqual(jones.colored_jones_polynomial([1]*3+[-2]*3, strands=3, color=N), multiply(a, b))

    def test_projector_graph_identities(self):
        A = Fraction(3, 2); delta = -A**2-A**-2
        for n in range(1, 6):
            p = colored.jones_wenzl_projector(n, A)
            self.assertEqual(diagram_product(p, p, delta), p)
            self.assertEqual(sum(c*delta**graph.trace_loops(idx) for idx, c in p.items()), colored.quantum_integer(n+1, A))
            for i in range(1, n):
                e = {graph.cup_cap(n, i): 1}
                self.assertFalse(diagram_product(e, p, delta))
                self.assertFalse(diagram_product(p, e, delta))

    def test_native_polynomial_arithmetic(self):
        P = jones._tl_tensor.PolyWrapper
        def wrap(p):
            low, high = min(p), max(p)
            return P(low, [p.get(e, 0) for e in range(low, high+1)])
        def terms(p): return {e: c for c, e in p.terms}
        rng = random.Random(202611041)
        for _ in range(100):
            a, b = [{e: rng.randrange(-10**20, 10**20) for e in rng.sample(range(-60, 61), 25)} for _ in range(2)]
            x, y = wrap(a), wrap(b)
            self.assertEqual(terms(x*y), multiply(a, b))
            self.assertEqual(terms(x+y), {e: a.get(e, 0)+b.get(e, 0) for e in a.keys() | b.keys() if a.get(e, 0)+b.get(e, 0)})
            self.assertFalse(x-x)
        for e in (-7, 0, 11):
            for power in (-10, -1, 0, 1, 20):
                self.assertEqual(terms(P(e, [-1])**power), {e*power: (-1)**(power % 2)})

    def test_compatibility_and_parameter_arithmetic(self):
        for A in (1, -1, 1j, -1j, Fraction(3, 2)):
            expected = colored.evaluate_polynomial(jones.colored_jones_polynomial([1]*3, strands=2, color=4), A)
            self.assertEqual(colored.colored_jones_numeric([1]*3, 2, 4, A), expected)
            self.assertEqual(colored.colored_jones_cabling([1]*3, 2, 4, A), expected)

    def test_native_import_independence(self):
        code = 'import sys; import rivin_jones as j; j.colored_jones_polynomial([1]*3,strands=2,color=5); assert not any(k=="sage" or k.startswith("sage.") for k in sys.modules); assert "tl_tensor" not in sys.modules'
        subprocess.run([sys.executable, '-B', '-c', code], check=True)


class Rejections(unittest.TestCase):
    def test_resource_limits(self):
        for option in ({'max_rows': 1}, {'max_transitions': 1}, {'max_color': 3}):
            with self.assertRaises(colored.ColoredResourceLimitError):
                jones.colored_jones_evaluate([1, -2]*2, strands=3, color=4, A=Fraction(2), **option)
        with self.assertRaises(colored.ColoredResourceLimitError):
            jones.colored_jones_polynomial([1]*3, strands=2, color=4, algorithm='r-matrix', max_polynomial_span=1)
        with self.assertRaises(colored.ColoredResourceLimitError):
            jones.colored_jones_evaluate([-1]*3, strands=2, color=4, A=Fraction(2), max_states=1)

    def test_bad_input(self):
        for options in [dict(strands=2, color=True), dict(strands=2, color=0), dict(strands=True, color=3),
                        dict(strands=3, color=3), dict(strands=2, color=3, normalize=1)]:
            with self.assertRaises((TypeError, ValueError)):
                jones.colored_jones_polynomial([1], **options)
        for A in (0, True, float('nan'), complex(float('inf'), 0)):
            with self.assertRaises((TypeError, ValueError)):
                jones.colored_jones_evaluate([1], strands=2, color=3, A=A)
        for word in ([0], [2], [True], [1.0]):
            with self.assertRaises((TypeError, ValueError)):
                jones.jones_polynomial(word, strands=2)

    def test_native_overflow_and_nonunits(self):
        P = jones._tl_tensor.PolyWrapper
        for call in [lambda: P(sys.maxsize, [1, 1]), lambda: P(sys.maxsize, [1])**2,
                     lambda: P(sys.maxsize, [1])*P(1, [1]), lambda: P(-1000001, [1])+P(0, [1]),
                     lambda: P(0, [2])**-1, lambda: P(0, [])**-1]:
            with self.assertRaises(ValueError): call()


if __name__ == '__main__': unittest.main()
