#!/usr/bin/env python3
"""Finite-color Arb evaluation; requires Sage for fields, not Jones computation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--word', required=True, help='JSON braid word')
    p.add_argument('--strands', required=True, type=int)
    p.add_argument('--colors', default='10,20,50')
    p.add_argument('--bits', type=int, default=512)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists(): p.error('choose a new output filename')
    if a.bits < 64: p.error('--bits must be at least 64')
    colors = [int(n) for n in a.colors.split(',')]
    if any(n < 2 for n in colors): p.error('colors must be at least 2')
    word = json.loads(a.word)
    report = dict(status='running', word=word, strands=a.strands, bits=a.bits, records=[],
                  convention='N-dimensional reduced colored Jones; q=A^-4, A=exp(-pi*i/(2N))',
                  scope='finite-N numerical enclosures, not a bound on convergence to the limiting volume',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with tempfile.TemporaryDirectory(prefix='rivin-jones-volume-') as folder:
        os.environ.update(DOT_SAGE=folder+'/sage', MPLCONFIGDIR=folder+'/mpl', XDG_CACHE_HOME=folder+'/cache')
        from sage.all import RealBallField, ComplexBallField
        from rivin_jones import colored_jones_evaluate
        R, C = RealBallField(a.bits), ComplexBallField(a.bits)
        for N in colors:
            start = time.perf_counter()
            A = C(0, -R.pi()/(2*N)).exp()
            value, stats = colored_jones_evaluate(word, strands=a.strands, color=N, A=A, return_info=True)
            if not abs(value) > 0: raise ArithmeticError('value enclosure contains zero; increase --bits or test exact vanishing')
            growth = 2*R.pi()/N*abs(value).log()
            report['records'].append(dict(color=N, value=str(value), growth=str(growth),
                lower=str(growth.lower().exact_rational()), upper=str(growth.upper().exact_rational()),
                seconds=time.perf_counter()-start, contraction=stats))
    report['status'] = 'passed'
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__': main()
