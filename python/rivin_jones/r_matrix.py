"""Sparse finite-dimensional U_v(sl2) braid contraction, v=A^2.

Basis vector i has weight h_i=n-2i, n=N-1. Positive crossings send
(i,j) to (j+k,i-k); negative crossings send (i,j) to (j-k,i+k).
All local coefficients are Laurent polynomials, computed without division.
Fix the first strand's input/output to its highest-weight vector and trace
the other strands with K. This reduced open-strand trace is regular at roots
of unity. Our framing factor is A^(-w*n*(n+2)). See the candidate derivation
and Huynh--Le, On the colored Jones polynomial and the Kashaev invariant, 2.1.
"""
from itertools import product

from .colored import ColoredResourceLimitError


class _Weights:
    def __init__(self, A, color, cache_limit):
        self.A = A
        self.n = color-1
        self.one = A**0
        self.zero = self.one-self.one
        self.cache_limit = cache_limit
        self.powers = {0: self.one}
        self.weights = {}
        # Symmetric Gaussian binomials: B(a,b)=v^-b B(a-1,b)+v^(a-b) B(a-1,b-1).
        self.binomial = [[self.one]]
        for a in range(1, color):
            previous = self.binomial[-1]
            self.binomial.append([self.one]+[
                self.power(-2*b)*previous[b]+self.power(2*(a-b))*previous[b-1]
                for b in range(1, a)]+[self.one])
        # Products of v^(m+r)-v^(-m-r), without quantum-factorial division.
        self.factors = []
        for m in range(color):
            row = [self.one]
            for r in range(1, color-m):
                row.append(row[-1]*(self.power(2*(m+r))-self.power(-2*(m+r))))
            self.factors.append(row)

    def power(self, exponent):
        value = self.powers.get(exponent)
        if value is None:
            value = self.A**exponent
            if len(self.powers) < self.cache_limit:
                self.powers[exponent] = value
        return value

    def coefficient(self, sign, i, j, k):
        key = (sign, i, j, k)
        value = self.weights.get(key)
        if value is not None:
            return value
        n = self.n
        if sign > 0:
            exponent = (n-2*i+2*k)*(n-2*j-2*k)+k*(k-1)
            value = self.power(exponent)*self.binomial[j+k][k]*self.factors[n-i][k]
        else:
            exponent = -(n-2*i)*(n-2*j)-k*(k-1)
            value = self.power(exponent)*self.binomial[i+k][k]*self.factors[n-j][k]
            if k % 2:
                value = -value
        if len(self.weights) < self.cache_limit:
            self.weights[key] = value
        return value


def _domains(word, n, target):
    """Backward interval over-approximation of feasible colors, using R support.

    This prunes impossible transitions, never coefficients by magnitude.
    Correlations between strand colors are deliberately not assumed.
    """
    low, high = list(target), list(target)
    after = [None]*len(word)
    for t in range(len(word)-1, -1, -1):
        g = word[t]
        p = abs(g)-1
        u, U, v, V = low[p], high[p], low[p+1], high[p+1]
        after[t] = (u, U, v, V)
        if g > 0:
            low[p], high[p] = v, min(n, U+V)
            low[p+1], high[p+1] = max(0, u+v-n), U
        else:
            low[p], high[p] = max(0, u+v-n), V
            low[p+1], high[p+1] = u, min(n, U+V)
    possible = all(a <= x <= b for a, x, b in zip(low, target, high))
    return after, possible


def _contract(word, strands, color, A, *, normalize=True, max_states=10000,
              max_rows=1000000, max_transitions=2000000, cache_limit=20000,
              prune=True):
    """Contract validated input; return a scalar in A's ring and work counters."""
    rows = color**(strands-1)
    if rows > max_rows:
        raise ColoredResourceLimitError(f'open-strand trace needs {rows} rows; limit is {max_rows}')
    weights = _Weights(A, color, cache_limit)
    n, one, zero = color-1, weights.one, weights.zero
    total = zero
    transitions = active = 0
    peak = 1
    for rest in product(range(color), repeat=strands-1):
        target = (0,)+rest
        if prune:
            domains, possible = _domains(word, n, target)
            if not possible:
                continue
        else:
            domains = [(0, n, 0, n)]*len(word)
        active += 1
        states = {target: one}
        for g, (u, U, v, V) in zip(word, domains):
            p = abs(g)-1
            sign = 1 if g > 0 else -1
            updated = {}
            for state, amplitude in states.items():
                i, j = state[p:p+2]
                if sign > 0:
                    lower = max(0, u-j, i-V)
                    upper = min(i, n-j, U-j, i-v)
                else:
                    lower = max(0, j-U, v-i)
                    upper = min(j, n-i, j-u, V-i)
                for k in range(lower, upper+1):
                    transitions += 1
                    if transitions > max_transitions:
                        raise ColoredResourceLimitError('direct colored transition limit exceeded')
                    coefficient = weights.coefficient(sign, i, j, k)
                    if coefficient == zero:
                        continue
                    pair = (j+k, i-k) if sign > 0 else (j-k, i+k)
                    output = state[:p]+pair+state[p+2:]
                    value = amplitude*coefficient
                    if output in updated:
                        value = updated[output]+value
                    elif len(updated) >= max_states:
                        raise ColoredResourceLimitError('direct colored frontier state limit exceeded')
                    if value == zero:
                        updated.pop(output, None)
                    else:
                        updated[output] = value
            states = updated
            peak = max(peak, len(states))
            if not states:
                break
        amplitude = states.get(target, zero)
        if amplitude != zero:
            trace_exponent = 2*((strands-1)*n-2*sum(rest))
            total = total+amplitude*weights.power(trace_exponent)
    writhe = sum(1 if g > 0 else -1 for g in word)
    if normalize:
        total = total*weights.power(-writhe*n*(n+2))
    elif (writhe*n) % 2:
        # Match the TL blackboard framing convention theta_n=(-1)^n A^(n(n+2)).
        total = -total
    stats = dict(algorithm='r-matrix', color=color, strands=strands, crossings=len(word),
                 trace_rows=rows, active_rows=active, transitions=transitions, peak_states=peak,
                 cached_weights=len(weights.weights), cached_powers=len(weights.powers),
                 backward_pruning=prune)
    return total, stats
