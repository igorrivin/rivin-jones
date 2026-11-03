# Planned reduced colored contraction

`plan_colored_jones(word, strands=s, color=N)` returns a frozen
`ColoredJonesPlan`. Its `evaluate(A)` preserves the scalar field;
`evaluate_many(As)` reuses the plan without storing answers. The optional
`backend='numpy'` batch path uses complex128 arrays. This is an explicit
method: automatic evaluation still uses the previous sparse braid sweep.

## Network and reduction

Each crossing uses the existing division-free `_Weights` coefficients. Wire
labels connect each output to the next input on that strand. Closing the
braid identifies final and initial labels. Strand zero is cut and its two
endpoints are fixed to state zero. Each other strand contributes its quantum
trace weight `A**(2*(N-1-2*i))` at its first input crossing. These factors are
inserted once, before any contraction, and the usual framing correction is
applied after the scalar sum. There is no division by a quantum dimension.

This is precisely the finite sum in [MATHEMATICS.md](MATHEMATICS.md), regrouped
according to the network. An index repeated within one crossing is traced
inside that leaf: only entries with equal repeated states contribute, and
that index is summed there. Every remaining wire occurs in exactly two
leaves. Input validation rejects multi-component closures before planning.

Scalar leaves enumerate only admissible `(i,j,k)` transitions. The relation
`i+j=i'+j'` is preserved, and sparse intermediate dictionaries store only
nonzero entries. Binary joins group entries by shared wire states and sum
matching products. This preserves rational, exact algebraic, native Laurent
polynomial and Arb arithmetic. Arb values are not discarded by a magnitude
threshold; only exact zero equality permits omission.

NumPy leaves use the same entry generator and coefficient formulas. Their
first axis indexes independent parameter values. Each binary contraction
permutes the other axes, reshapes to matrices, and performs batched matrix
multiplication. The batch index is never contracted. Output storage is always
complex128, so integer/real input arrays cannot truncate complex answers.
This path produces numerical approximations; it does not provide Arb-style
enclosures or an error certificate.

## Planning and reuse

The default is Cotengra's greedy contraction path. `optimize` can also be an
unsliced Cotengra optimizer or a standard binary linear path. `'sequential'`
uses a deterministic left fold for comparisons. A supplied hyperoptimizer
controls its own search-time/repeat/parallelism budget. Pre-sliced optimizer
trees are rejected; slicing is performed by the checked plan builder.

Only immutable tuples of topology, path steps and sliced wire labels are
retained. A plan is specific to its braid, width and color. Its coefficients
and bounded caches are reconstructed per evaluation or NumPy chunk. Changing
A or normalization is supported; changing color requires a new plan. In
particular, a Volume Conjecture sequence changes the network dimensions.

## Resource controls

Planning defaults:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `max_color` | 256 | Representation dimension |
| `max_crossings` | 1,000 | Input crossing count |
| `max_entries` | 1,000,000 | Estimated live scalar slots |
| `max_slices` | 4,096 | Maximum slice assignments |
| `max_operations` | 100,000,000 | Conservative scalar-update work estimate per point |
| `cache_limit` | 20,000 | Maximum entries in each coefficient/power cache |

The live-slot estimate sums live dense tensor sizes, an output allocation,
possible contiguous reshape copies, and coefficient/cache storage. It is
deliberately conservative for the sparse scalar backend. Cache reserves are
bounded both by `cache_limit` and the number of possible coefficients/powers.
Setting `cache_limit=0` can help when requesting a very small tensor budget.

When the estimate exceeds `max_entries`, the planner fixes a wire index,
removes it from the live tensor shapes, and recomputes the bound. It chooses
the best resulting memory bound, then work estimate, with a deterministic
tie-break. For k fixed wires it evaluates all `N**k` assignments and adds
their results. This is an exact decomposition, without sampled slices,
truncated states or altered normalization. Unsatisfiable budgets raise
`ColoredResourceLimitError` before coefficients are constructed.

The NumPy batch API defaults to `batch_size=128`, `max_bytes=64 MiB`, and
`max_points=100000`. It limits input consumption, validates every parameter,
and shrinks chunks according to the plan's slot bound, complex128 item size,
and input/output arrays. A budget that cannot hold one point raises instead
of attempting an oversized allocation. Floating overflow or nonfinite
parameters/results fail explicitly.

These are workload and array-payload estimates, **not process RSS caps**.
Python containers, BLAS workspaces and the variable bit sizes of exact field
elements are outside that model. Slicing can increase total work sharply;
the slice and work limits prevent that increase from being silently accepted.
Dense large-color networks may be unsuitable even with slicing. The existing
sparse sweep and its backward pruning remain available for those cases.

## Validation and performance scope

Portable tests replay all 176 frozen colored coefficient fixtures through the
planned contraction using native exact Laurent arithmetic. Other checks cover
Markov stabilization, framing, alternative paths, forced slicing, repeated
evaluation, invalid inputs, preallocation guards and integer/complex batches
crossing the 128-point boundary. Optional Sage tests compare fresh colored
coefficients at exact roots and check sliced/unsliced Arb results against the
independent figure-eight formula.

The CPU study uses six fixed knot/color cases, compares the previous sweep
with scalar and NumPy plan execution, and records planning separately. Arb
comparisons use 256-bit enclosures. Scalar timings are medians of five
post-warm-up calls; NumPy 32-point batches use three calls. The batch loop and
polynomial-first baselines are single measurements. Narrow cases favor the
existing engine; wide cases can benefit substantially. These data do not
establish a universal speedup, a Sage comparison, or GPU performance.

No GPU backend is implemented in this change. CPU correctness, root behavior,
memory accounting and batching are established first. GPU kernels, precision
policy, transfer/synchronization costs and device memory require separate
validation before being offered as a supported backend.
