"""Independent coefficient, topology, slicing, and numeric batch regressions."""
from dataclasses import FrozenInstanceError
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import unittest

import rivin_jones as jones
from rivin_jones.colored import ColoredResourceLimitError, evaluate_polynomial


class PlannedContraction(unittest.TestCase):
    def test_all_frozen_colored_coefficients(self):
        cases = json.loads((Path(__file__).parent/'data/colored.json').read_text())
        for c in cases:
            with self.subTest(name=c['name'], color=c['color'], word=c['word']):
                plan = jones.plan_colored_jones(c['word'], strands=c['strands'], color=c['color'])
                poly = plan.evaluate(jones._tl_tensor.PolyWrapper(1, [1]))
                self.assertEqual({e: v for v, e in poly.terms},
                                 {int(e): v for e, v in c['coefficients'].items()})

    def test_reuse_slices_order_and_framing(self):
        word = [1, -2]*2
        params = [Fraction(3, 2), Fraction(4, 3), 1j, -1]
        coefficients = jones.colored_jones_polynomial(word, strands=3, color=3)
        for optimize in ('greedy', 'sequential'):
            plan = jones.plan_colored_jones(word, strands=3, color=3,
                                           optimize=optimize, max_entries=100, cache_limit=0)
            self.assertGreater(plan.slices, 1)
            self.assertLessEqual(plan.peak_entries, 100)
            for A, actual in zip(params, plan.evaluate_many(params)):
                self.assertEqual(actual, evaluate_polynomial(coefficients, A))
        self.assertRaises(FrozenInstanceError, setattr, plan, 'color', 4)
        for word, strands in [([1]*3, 2), ([-1]*3, 2), ([1]*3+[2], 3), ([1]*3+[-2], 3)]:
            for N in (2, 3, 4):
                plan = jones.plan_colored_jones(word, strands=strands, color=N)
                for normalize in (True, False):
                    expected = jones.colored_jones_evaluate(word, strands=strands, color=N,
                        A=Fraction(3, 2), normalize=normalize)
                    self.assertEqual(plan.evaluate(Fraction(3, 2), normalize=normalize), expected)
                    self.assertEqual(jones.colored_jones_evaluate(word, strands=strands, color=N,
                        A=Fraction(3, 2), normalize=normalize, algorithm='tensor'), expected)

    def test_validation_and_preallocation_guards(self):
        for word, strands, color in [([0], 2, 2), ([True], 2, 2), ([2], 2, 2),
                                      ([1], 1, 2), ([1], 3, 2), ([1], 2, True)]:
            with self.assertRaises((TypeError, ValueError)):
                jones.plan_colored_jones(word, strands=strands, color=color)
        for options in ({'max_entries': 1}, {'max_operations': 1}, {'max_color': 2},
                        {'max_crossings': 3}, {'max_entries': 100, 'max_slices': 1, 'cache_limit': 0}):
            with self.assertRaises(ColoredResourceLimitError):
                jones.plan_colored_jones([1, -2]*2, strands=3, color=3, **options)
        plan = jones.plan_colored_jones([1]*3, strands=2, color=3)
        for A in (True, 0, float('nan'), complex(float('inf'), 0)):
            with self.assertRaises((ValueError, TypeError)):
                plan.evaluate(A)
        for kwargs in ({'batch_size': 0}, {'batch_size': True}, {'max_bytes': 0},
                       {'max_points': 0}, {'normalize': 1}, {'backend': 'unknown'}):
            with self.assertRaises((ValueError, TypeError)):
                plan.evaluate_many([1], **kwargs)
        with self.assertRaises(ColoredResourceLimitError):
            plan.evaluate_many(iter([1, 2, 3]), max_points=2)
        self.assertEqual(plan.evaluate_many([]), [])
        self.assertEqual(jones.plan_colored_jones([], strands=1, color=200,
                                                max_entries=1).evaluate(1j), 1)


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'optional NumPy tensor extra')
class NumpyContraction(unittest.TestCase):
    def test_batches_preserve_fractional_and_complex_results(self):
        import numpy as np
        plan = jones.plan_colored_jones([-1]*3, strands=2, color=3)
        for points in (np.full(129, 2, dtype=int), np.full(129, 0.9+0.2j),
                       np.array([1, -1, 1j, -1j, 1.1])):
            expected = np.array([complex(x) for x in plan.evaluate_many(points)])
            budget = 32*len(points)+16*plan.peak_entries*2
            actual, info = plan.evaluate_many(points, backend='numpy', max_bytes=budget, return_info=True)
            self.assertEqual(actual.dtype, np.complex128)
            self.assertTrue(np.allclose(actual, expected, rtol=1e-10, atol=1e-10))
            self.assertLessEqual(info['batch_size'], 2)
            self.assertLessEqual(info['estimated_peak_bytes'], budget)
        self.assertEqual(plan.evaluate_many([], backend='numpy').shape, (0,))
        with self.assertRaises(ColoredResourceLimitError):
            plan.evaluate_many([1], backend='numpy', max_bytes=1)
        with self.assertRaises((ValueError, TypeError)):
            plan.evaluate_many([[1, 2]], backend='numpy')
        with self.assertRaises(TypeError):
            plan.evaluate_many(np.array([True]), backend='numpy')
        with self.assertRaises(TypeError):
            jones.plan_colored_jones([1], strands=2, color=np.bool_(True))
        for A in (np.float32('nan'), np.complex64(complex(float('inf'), 0))):
            with self.assertRaises(ValueError):
                plan.evaluate(A)

    def test_slicing_and_roots(self):
        import numpy as np
        for N in (2, 3, 4):
            plan = jones.plan_colored_jones([1, -2]*2, strands=3, color=N,
                max_entries=150, cache_limit=0, max_operations=100000000)
            points = np.array([np.exp(-1j*np.pi/(2*N)), 0.9+0.1j, 1.1-0.2j])
            expected = plan.evaluate_many(points)
            actual = plan.evaluate_many(points, backend='numpy', batch_size=2)
            self.assertTrue(np.allclose(actual, expected, rtol=1e-9, atol=1e-9))


if __name__ == '__main__':
    unittest.main()
