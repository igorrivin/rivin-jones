# Direct colored contraction: conventions and reduction

This implementation uses representation dimension N, n=N-1, and
`v=A^2`. Its Jones variable is `q=t=A^-4=v^-2`. The reduced unknot is 1.
Input is an explicitly sized braid whose closure has one component.

The construction is the finite-dimensional quantum-group braid action followed
by a partial quantum trace. See [Huynh–Lê, Section 2.1.4](https://letu.math.gatech.edu/Papers/Jones_Kashaev.pdf)
for the highest-weight partial-trace construction and
[Park, Section 2](https://www.categorical.center/Large_color_R_matrix_for_knot_complements_and_strange_identities.pdf)
for colored R matrices and the open-strand trace. Parameter, weight ordering
and crossing conventions differ between references. The formulas below specify
the convention used in this candidate completely; the independent cable
calculation fixes its identification with the workspace's Jones polynomial.

## Local matrices

Basis vector i, for 0 <= i <= n, has weight h_i=n-2i. The symmetric Gaussian
binomial B(a,b) is a Laurent polynomial in v, computed by

\[
B(a,b)=v^{-b}B(a-1,b)+v^{a-b}B(a-1,b-1),\qquad B(a,0)=B(a,a)=1.
\]

The positive crossing takes (i,j) to (j+k,i-k), for
0 <= k <= min(i,n-j), with coefficient

\[
A^{(n-2i+2k)(n-2j-2k)+k(k-1)}
B(j+k,k)\prod_{r=1}^{k}
\left(A^{2(n-i+r)}-A^{-2(n-i+r)}\right).
\]

Its inverse takes (i,j) to (j-k,i+k), for
0 <= k <= min(j,n-i), with coefficient

\[
(-1)^k A^{-(n-2i)(n-2j)-k(k-1)}
B(i+k,k)\prod_{r=1}^{k}
\left(A^{2(n-j+r)}-A^{-2(n-j+r)}\right).
\]

These expressions contain no divisions by quantum integers or factorials.
They therefore specialize at any nonzero A, even when the root order differs
from the representation dimension. Algebraic cancellations must not be
replaced by a threshold on coefficient magnitude.

## Reduced closure and framing

For an s-strand braid w, write R(w) for its action on the tensor power. Fix
strand zero's input and output to its highest-weight vector i_0=0. Trace the
remaining strands with K, which has eigenvalue A^(2h_i). The reduced invariant is

\[
J_N(\widehat w;A)=A^{-\operatorname{wr}(w)n(n+2)}
\sum_{i_1,\ldots,i_{s-1}=0}^{n}
A^{2\sum_{r=1}^{s-1}(n-2i_r)}
R(w)_{(0,i_1,\ldots,i_{s-1}),(0,i_1,\ldots,i_{s-1})}.
\]

At generic A, the partially traced knot endomorphism is scalar; its value on
the highest-weight vector gives the invariant divided by the quantum dimension.
The displayed finite Laurent expression extends that reduced invariant to
roots of unity. It computes no quotient of two specialized zero values.

The ambient correction above has **no additional writhe sign**. For a knot
closure, writhe and s-1 have the same parity; the signs introduced by translating
the usual quantum trace to the TL circle convention cancel in the normalized
invariant. To return the workspace's blackboard-framed reduced value, multiply
the uncorrected trace by (-1)^(n*wr(w)). Equivalently,

\[
J_N^{\rm framed}=(-1)^{n\operatorname{wr}(w)}
A^{\operatorname{wr}(w)n(n+2)} J_N.
\]

Color 2 is the ordinary reduced Jones polynomial. A normalized knot result
must have only A-exponents divisible by four. Neither the trace reduction nor
the sign translation here defines a multi-component colored-link convention.

## Contraction and bounds

Each diagonal trace row starts as one sparse tensor-product basis state.
Crossings change only the two adjacent strand colors. Equal resulting states
are combined exactly. Local coefficients and powers have bounded per-call
caches; no answers are cached between calls.

Backward propagation from the target tuple gives interval over-approximations
of feasible colors. For output ranges [u,U] and [v,V], a positive crossing has
input ranges [v,min(n,U+V)] and [max(0,u+v-n),U]. A negative crossing has input
ranges [max(0,u+v-n),V] and [u,min(n,U+V)]. These follow directly from the
transition formulas and conserve i+j. Ignoring correlations enlarges the
allowed set; it cannot remove a genuine path. Forward transitions outside
these ranges cannot contribute to the requested diagonal entry.

For exact polynomials, the scalar arithmetic uses the existing native GMP
Laurent representation. Each local coefficient has exponents within
[-6n^2,6n^2], a conservative bound obtained by bounding the displayed monomial,
Gaussian binomial and product separately. The trace adds at most 2(s-1)n in
either direction. Thus `12*c*n^2+4*(s-1)*n` bounds the support span before
framing translation. The caller's span guard is checked before contraction;
the native arithmetic bridge also checks spans and exponent overflow.

Trace rows grow as N^(s-1); states and transitions can also grow sharply.
The figure-eight presentation happens to prune to one live state per row,
but this is not a general complexity bound. Runtime guards reject an overrun
without returning a partial invariant.

## Verification boundaries

The audit checks the local inverse and Yang–Baxter equations on every basis
vector in small representations, over both rational and cyclotomic fields.
It compares complete polynomials with cable/Sage references, tests both Markov
signs, mirror and framing behavior, compares pruned and unpruned contractions,
and checks exact specialized values at several independent root orders.

For the Volume Conjecture, use A_N=exp(-pi*i/(2N)). Complex ball arithmetic
encloses the finite-N value. Increasing precision and comparison with the
independent figure-eight product formula test numerical evaluation separately
from polynomial equality. Those intervals do not bound the error to the
N-to-infinity limit, and these finite checks are not a universal implementation
proof or a guarantee of performance for all braids.
