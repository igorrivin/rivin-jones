"""Optional live Sage and Arb checks, enabled explicitly in a Sage environment."""
import os
import tempfile
import unittest

import rivin_jones as jones
from rivin_jones import colored


@unittest.skipUnless(os.environ.get('RIVIN_JONES_TEST_SAGE') == '1', 'set RIVIN_JONES_TEST_SAGE=1 under Sage Python')
class SageOracles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache = tempfile.TemporaryDirectory(prefix='rivin-jones-sage-')
        cls.old = {k: os.environ.get(k) for k in ('DOT_SAGE', 'MPLCONFIGDIR', 'XDG_CACHE_HOME')}
        os.environ.update(DOT_SAGE=cls.cache.name+'/sage', MPLCONFIGDIR=cls.cache.name+'/mpl', XDG_CACHE_HOME=cls.cache.name+'/cache')

    @classmethod
    def tearDownClass(cls):
        for key, value in cls.old.items():
            if value is None: os.environ.pop(key, None)
            else: os.environ[key] = value
        cls.cache.cleanup()

    def test_fresh_sage_polynomials_and_exact_roots(self):
        from sage.all import BraidGroup, CyclotomicField
        for word, strands in [([1]*3, 2), ([-1]*3, 2), ([1, -2]*2, 3), ([1]*3+[-2], 3)]:
            for N in (2, 3, 4, 5):
                expected = {-4*int(e): int(c) for e, c in BraidGroup(strands)(list(word)).colored_jones_polynomial(N).dict().items() if c}
                self.assertEqual(jones.colored_jones_polynomial(word, strands=strands, color=N), expected)
                plan = jones.plan_colored_jones(word, strands=strands, color=N)
                for order in sorted({4*N, 8, 12}):
                    A = CyclotomicField(order).gen()**-1
                    self.assertEqual(jones.colored_jones_evaluate(word, strands=strands, color=N, A=A), colored.evaluate_polynomial(expected, A))
                    self.assertEqual(plan.evaluate(A), colored.evaluate_polynomial(expected, A))

    def test_planned_arb_and_sliced_roots(self):
        from sage.all import RealBallField, ComplexBallField
        R, C = RealBallField(512), ComplexBallField(512)
        for N, budget in ((3, 100), (6, 1000), (10, 1000000)):
            plan = jones.plan_colored_jones([1, -2]*2, strands=3, color=N,
                max_entries=budget, cache_limit=0)
            A = C(0, -R.pi()/(2*N)).exp()
            value = plan.evaluate(A)
            term = total = R(1)
            for k in range(1, N):
                term *= (2*(R.pi()*k/N).sin())**2
                total += term
            self.assertTrue(value.real().overlaps(total))
            self.assertTrue(value.imag().contains_zero())
            self.assertGreater(value.real().accuracy(), 450)

    def test_arb_roots_against_independent_formula(self):
        from sage.all import RealBallField, ComplexBallField
        for N in (10, 50, 200):
            bits = max(256, 4*N+64)
            R, C = RealBallField(bits), ComplexBallField(bits)
            A = C(0, -R.pi()/(2*N)).exp()
            value = jones.colored_jones_evaluate([1, -2]*2, strands=3, color=N, A=A)
            term = total = R(1)
            for k in range(1, N):
                term *= (2*(R.pi()*k/N).sin())**2
                total += term
            self.assertTrue(value.real().overlaps(total))
            self.assertTrue(value.imag().contains_zero())
            self.assertTrue(abs(value) > 0)


if __name__ == '__main__': unittest.main()
