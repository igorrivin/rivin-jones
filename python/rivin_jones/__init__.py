from . import _tl_tensor

import colorsys
import dataclasses
import math
import functools
import operator
import cotengra


def colored_jones_polynomial(word, *, strands, color, **options) -> dict[int, int]:
    """Reduced N-dimensional colored Jones of a knot; q=A^-4. See .colored."""
    from .colored import colored_jones_polynomial as compute
    return compute(word, strands=strands, color=color, **options)


def colored_jones_evaluate(word, *, strands, color, A, **options):
    """Reduced colored evaluation by a direct open-strand trace, including at roots."""
    from .colored import colored_jones_evaluate as compute
    return compute(word, strands=strands, color=color, A=A, **options)


def plan_colored_jones(word, *, strands, color, **options):
    """Build a reusable reduced colored tensor plan; see .tensor for limits."""
    from .tensor import plan_colored_jones as compile_plan
    return compile_plan(word, strands=strands, color=color, **options)


def jones_polynomial(word, *, strands, method="auto") -> dict[int, int]:
    """Exact reduced Jones coefficients in A, with t=A^-4 and explicit width.

    auto uses the packed native braid frontier through eight strands. Through
    twelve strands it compares a frontier state bound with fresh greedy work;
    larger widths use greedy contraction. Frontier arithmetic uses i128
    only under a proved coefficient bound, otherwise GMP. Nothing is cached
    across calls and no external Jones engine is used.
    """
    if isinstance(strands, bool):
        raise TypeError("strands must be an integer, not bool")
    strands = operator.index(strands)
    if not 1 <= strands <= 128:
        raise ValueError("strands must be between 1 and 128")
    generators = []
    for g in word:
        if isinstance(g, bool):
            raise TypeError("braid generators must be integers, not bool")
        g = operator.index(g)
        if not 1 <= abs(g) < strands:
            raise ValueError("braid generators must satisfy 1 <= abs(g) < strands")
        generators.append(g)
    if method == "auto" and not generators:
        # Preserve every declared idle component in the reduced unlink value.
        power = strands-1
        sign = -1 if power % 2 else 1
        return {2*power-4*j: sign*math.comb(power, j) for j in range(power+1)}
    network = None
    path = None
    if method == "auto":
        method = "frontier" if strands <= 8 else "greedy"
        if 8 < strands <= 12:
            network = TLTensorNetwork.from_word(generators, strands=strands)
            info = network.contract_info("greedy")
            # Missing generators separate independent TL subalgebras. Their
            # Catalan dimensions give a bound on the reachable frontier states.
            used = {abs(g) for g in generators}
            cuts = [0] + [i for i in range(1, strands) if i not in used] + [strands]
            states = math.prod(_catalan(b-a) for a, b in zip(cuts, cuts[1:]))
            frontier_work = 4*len(generators)*states
            if frontier_work < 10**info.cost:
                method = "frontier"
            else:
                path = info.path
    if method in ("frontier", "frontier-gmp"):
        values, _ = _tl_tensor.braid_jones(generators, strands,
                                          "gmp" if method == "frontier-gmp" else "auto")
        return {exponent: coefficient for coefficient, exponent in values}
    if method not in ("greedy", "sequential"):
        raise ValueError("method must be auto, frontier, frontier-gmp, greedy or sequential")
    if network is None:
        network = TLTensorNetwork.from_word(generators, strands=strands)
    answer = network.contract(optimize=path if path is not None else method, parallel=False)
    if len(answer.tensors) != 1 or answer.tensors[0].bonds:
        raise RuntimeError("contraction did not produce a scalar")
    terms = answer.tensors[0].terms
    if len(terms) != 1 or terms[0][0].idx or terms[0][0].loops != 1:
        raise RuntimeError("unexpected reduced scalar diagram")
    return {exponent: coefficient for coefficient, exponent in terms[0][1]}

def _ach(st):
    if not isinstance(st, str):
        st = str(st)
    hs = hash(st)
    h = (hs % 100) / 100
    s = 0.3 + 0.5 * (((hs // 100) % 100) / 100)
    v = 0.8 + 0.1 * (((hs // 10000) % 100) / 100)
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    chex = f"#{r:02x}{g:02x}{b:02x}"
    return f"<b style='color: {chex};'>{st}</b>"

LPoly = list[tuple[int, int]]

class TLDiagram:
    def __init__(self, idx: list[int], loops: int = 0):
        assert len(idx) <= 256 and len(idx) % 2 == 0
        assert all(0 <= i < len(idx) for i in idx)
        assert loops >= 0

        self._inner = _tl_tensor.TLDiagramInner(idx, loops)

    @staticmethod
    def eye(n: int, loops: int = 0) -> 'TLDiagram':
        assert n >= 0
        return TLDiagram(list(range(n, 2*n)) + list(range(n)), loops=loops)
    
    @staticmethod
    def gen(n: int, k: int, loops: int = 0) -> 'TLDiagram':
        assert n >= 0 and 0 <= k < n - 1
        idx = list(range(n, 2*n)) + list(range(n))
        idx[k] = k + 1
        idx[k + 1] = k
        idx[k + n] = k + n + 1
        idx[k + n + 1] = k + n
        return TLDiagram(idx, loops=loops)
    
    @staticmethod
    def _from_inner(inner) -> 'TLDiagram':
        res = TLDiagram.__new__(TLDiagram)
        res._inner = inner
        return res

    @property
    def idx(self):
        return list(self._inner.idx)
    
    @property
    def loops(self):
        return self._inner.loops
    
    @property
    def n(self):
        return len(self._inner.idx) // 2
    
    def tensor(self, other: 'TLDiagram') -> 'TLDiagram':
        inner = self._inner.tensor(other._inner)
        return TLDiagram._from_inner(inner)
    
    def __or__(self, other: 'TLDiagram') -> 'TLDiagram':
        return self.tensor(other)
    
    def trace(self, i: int, j: int):
        n = self.n
        assert 0 <= i <= 2*n and 0 <= j <= 2*n
        assert (i < n and j >= n) or (j < n and i >= n)
        self._inner.trace(i, j)

    def multiply(self, other: 'TLDiagram') -> 'TLDiagram':
        assert self.n == other.n

        inner = self._inner.tensor(other._inner)
        n = self.n
        for i in range(n):
            inner.trace(n, 2*n - i)
        
        return TLDiagram._from_inner(inner)
    
    def __matmul__(self, other: 'TLDiagram') -> 'TLDiagram':
        return self.multiply(other)
    
    def __mul__(self, other: LPoly) -> 'TLTensor':
        if isinstance(other, list):
            return TLTensor([(self, other)], list(range(2*self.n)))
        elif isinstance(other, int):
            return TLTensor([(self, [(other, 0)])], list(range(2*self.n)))
        else:
            raise ValueError()
    
    def __rmul__(self, other: LPoly) -> 'TLTensor':
        return self * other

    def draw(self, *, scale: float = 1.0, fig_kwargs: dict = None, ax_kwargs: dict = None):
        import matplotlib.pyplot
        import matplotlib.patches
        import matplotlib.lines
        import matplotlib.path

        maxc1 = max(((self.idx[i] - i) if self.idx[i] < self.n else 1 for i in range(self.n)), default=0)
        maxc2 = max(((self.idx[i] - i) if self.idx[i] >= self.n else 1 for i in range(self.n, 2*self.n)), default=0)
        max1 = max((abs((self.idx[i] % self.n) - i) for i in range(self.n)), default=0)
        max2 = max((abs((self.idx[i] % self.n) - (i % self.n)) for i in range(self.n, 2*self.n)), default=0)
        sep = max1/2 + max2/2 + 0.5
        boff = max(sep/2 - maxc1/3, sep/2 - maxc2/3)

        if fig_kwargs is None:
            fig_kwargs = {}
        if ax_kwargs is None:
            ax_kwargs = {}

        if 'dpi' not in fig_kwargs:
            fig_kwargs['dpi'] = 100
        if 'figsize' not in fig_kwargs:
            if self.loops > 0:
                fig_kwargs['figsize'] = [scale * (sep + 1) / 2, scale * (max(self.n, 2) - 0.8) / 2]
            else:
                fig_kwargs['figsize'] = [scale * sep / 2, scale * (max(self.n, 2) - 0.8) / 2]

        fig = matplotlib.pyplot.figure(**fig_kwargs)
        ax = fig.add_axes((0, 0, 1, 1), **ax_kwargs)

        if self.n >= 2:
            for i in range(self.n):
                if self.idx[i] < self.n and i < self.idx[i]:
                    ax.add_artist(matplotlib.patches.Arc([0.0, self.n - 1 - (i + self.idx[i]) / 2], (self.idx[i] - i), (self.idx[i] - i), theta1=270, theta2=90))
                elif self.idx[i] >= self.n:
                    ax.add_artist(matplotlib.patches.PathPatch(matplotlib.path.Path(
                        [(0, self.n - 1 - i), (maxc1/3, self.n - 1 - i), (maxc1/3+boff, self.n - 1 - i), (sep-maxc2/3-boff, 2*self.n - 1 - self.idx[i]), (sep-maxc2/3, 2*self.n - 1 - self.idx[i]), (sep, 2*self.n - 1 - self.idx[i])],
                        [matplotlib.path.Path.MOVETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4, matplotlib.path.Path.LINETO]
                    ), fc='none'))
            for i in range(self.n, 2*self.n):
                if self.idx[i] >= self.n and i < self.idx[i]:
                    ax.add_artist(matplotlib.patches.Arc([sep, 2*self.n - 1 - (i + self.idx[i]) / 2], (self.idx[i] - i), (self.idx[i] - i), theta1=90, theta2=270))
        elif self.n == 1:
            ax.add_artist(matplotlib.lines.Line2D([0, sep], [0.5, 0.5], color='k'))

        if self.loops > 0:
            ax.add_artist(matplotlib.patches.Circle([sep + 0.6, max(self.n, 2) / 2 - 1/2], 0.25, fc='none', ec='k'))
            ax.annotate(str(self.loops), (sep + 0.6, max(self.n, 2) / 2 - 1/2), ha='center', va='center', fontsize=10*scale)
            ax.set_xlim(0, sep + 1)
        else:
            ax.set_xlim(0, sep)
        ax.set_ylim(-0.1, max(self.n, 2) - 0.9)
        ax.set_axis_off()
        return fig

    def __str__(self):
        if self._inner.loops > 0:
            return f"TLDiagram({list(self._inner.idx)}, loops={self._inner.loops})"
        else:
            return f"TLDiagram({list(self._inner.idx)})"
    
    def __repr__(self):
        return str(self)
    
    def _ipython_display_(self):
        try:
            import matplotlib
        except:
            import IPython.display
            IPython.display.display_html(self._repr_html_())
            return
        self.draw()

    def _repr_html_(self):
        try:
            import matplotlib.pyplot
        except:
            if self._inner.loops > 0:
                return f"{_ach('TLDiagram')}([{', '.join(map(_ach, self._inner.idx))}], {_ach('loops')}={_ach(self._inner.loops)})"
            else:
                return f"{_ach('TLDiagram')}([{', '.join(map(_ach, self._inner.idx))}])"
        return _draw_diagram_to_html(tuple(self.idx), self.loops)
    
    def __eq__(self, other):
        return isinstance(other, TLDiagram) and self.loops == other.loops and self.idx == other.idx

@functools.lru_cache
def _draw_diagram_to_html(idx: tuple[int], loops: int):
    import io
    import base64
    import matplotlib.pyplot
    diag = TLDiagram(idx, loops)
    fig = diag.draw(scale=0.25)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0.05)
    matplotlib.pyplot.close(fig)
    img = base64.b64encode(buf.getvalue()).decode('utf-8')
    return f'<img src="data:image/png;base64,{img}" style="vertical-align: middle;"/>'

def _poly_to_rust(poly):
    start, end = min(power for _, power in poly), max(power for _, power in poly)
    coeffs = [0] * (end - start + 1)
    for coeff, power in poly:
        coeffs[power - start] = coeff
    return start, coeffs

class TLTensor:
    _BRACKET_LOOP_POLY = _tl_tensor.PolyWrapper(-2, [-1, 0, 0, 0, -1])

    def __init__(self, terms: list[tuple[TLDiagram, LPoly]], bonds: list[int]):
        assert all(b >= 0 for b in bonds)
        assert all(diag.n == terms[0][0].n for diag, _ in terms)
        assert 2*terms[0][0].n == len(bonds)

        self._inner = _tl_tensor.TLTensorInner([
            (diag._inner, *_poly_to_rust(poly))
            for diag, poly in terms
        ], bonds)

    @staticmethod
    def crossing(bonds: list[int], over=True):
        assert len(bonds) == 4
        if over:
            return TLTensor([(TLDiagram.eye(2), [(1, 1)]), (TLDiagram.gen(2, 0), [(1, -1)])], bonds)
        else:
            return TLTensor([(TLDiagram.eye(2), [(1, -1)]), (TLDiagram.gen(2, 0), [(1, 1)])], bonds)
        
    @staticmethod
    def scalar(poly: LPoly):
        return TLTensor([(TLDiagram([]), poly)], [])
        
    @staticmethod
    def _from_inner(inner) -> 'TLTensor':
        res = TLTensor.__new__(TLTensor)
        res._inner = inner
        return res
        
    @property
    def bonds(self) -> list[int]:
        return self._inner.bonds
    
    @property
    def terms(self) -> list[tuple[TLDiagram, LPoly]]:
        return [
            (TLDiagram._from_inner(diag), poly)
            for diag, poly in self._inner.terms
        ]
    
    @property
    def nterms(self) -> int:
        return self._inner.nterms
    
    def contract(self, 
        other: 'TLTensor', 
        poly: list[tuple[int, int]] = None, 
        progbar: bool = False,
        parallel: int = False
    ) -> 'TLTensor':
        cost = self.nterms * other.nterms
        prog = _tl_tensor.ProgressWrapper(cost, progbar)
        if parallel == "auto":
            pool = _tl_tensor.PoolWrapper(0)
        elif parallel == False:
            pool = _tl_tensor.PoolWrapper(-1)
        else:
            pool = _tl_tensor.PoolWrapper(parallel)
        if poly is None:
            poly = TLTensor._BRACKET_LOOP_POLY
        else:
            poly = _tl_tensor.PolyWrapper(*_poly_to_rust(poly))
        return self._contract(other, prog, pool, poly)
    
    def _contract(self, other: 'TLTensor', prog, pool, poly) -> 'TLTensor':
        return TLTensor._from_inner(self._inner.contract(other._inner, prog, pool, poly))
    
    def __xor__(self, other: 'TLTensor') -> 'TLTensor':
        return self.contract(other)

    def reindex(self, bonds: list[int]) -> 'TLTensor':
        return TLTensor(self.terms, bonds)

    def __add__(self, other: 'TLTensor') -> 'TLTensor':
        assert isinstance(other, TLTensor)
        assert other.bonds == self.bonds

        inner = self._inner.add(other._inner)
        return TLTensor._from_inner(inner)
    
    def __str__(self):
        return f"TLTensor({self.terms}, {self.bonds})"
    
    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        return isinstance(other, TLTensor) and self.bonds == other.bonds and self.terms == other.terms
    
    def _repr_html_(self):
        s = "<samp style='font-size: 12px;'>"
        s += "<details>"
        s += f"<summary>{_ach('TLTensor')}(bonds=[{', '.join(map(_ach, self.bonds))}], terms={_ach(len(self.terms))})</summary><table>"
        terms = self.terms
        for start in range(0, len(terms), 5):
            s += "<tr>"
            for idx, (d, p) in enumerate(terms[start:start+5], start=start):
                poly_str = ''
                for i, (c, pow) in enumerate(p):
                    if i > 0:
                        if c >= 0:
                            poly_str += ' + '
                        else:
                            poly_str += ' - '
                    elif c < 0:
                        poly_str += '-'
                    if pow == 0:
                        poly_str += f'{_ach(abs(c))}'
                    elif pow == 1:
                        poly_str += f'{_ach(abs(c))}x'
                    else:
                        poly_str += f"{_ach(abs(c))}x^{_ach(pow)}"
                if idx > 0:
                    s += "<td>+</td>"
                else:
                    s += "<td></td>"
                s += f"<td>{d._repr_html_()}&nbsp;&middot;&nbsp;({poly_str})</td>"
                if idx >= 99:
                    s += "</tr>"
                    s += "<p>...</p>"
                    break
            else:
                s += "</tr>"
                continue
            break
        s += "</table></details></samp>"
        return s
    
@dataclasses.dataclass
class ContractionInfo:
    cost: float
    max_intermediate: float
    path: list[tuple[int, int]]

    def _repr_html_(self):
        s = "<samp style='font-size: 12px;'>"
        s += "<details>"
        s += f"<summary>{_ach('ContractionInfo')}(log10(cost)={_ach(f'{self.cost:.3g}')}, log2(max_intermediate)={_ach(f'{self.max_intermediate:.3g}')})</summary>"
        s += 'path=[<br/>'
        for i in range(0, len(self.path), 10):
            s += "&nbsp;&nbsp;" + ", ".join(f"({_ach(i)}, {_ach(j)})" for i, j in self.path[i:i+10])
            if i < len(self.path) - 10:
                s += ',<br/>'
            else:
                s += '<br/>'
        s += "]</details></samp>"
        return s

class TLTensorNetwork:
    def __init__(self, tensors: list[TLTensor]):
        self.tensors = tensors

    def __str__(self):
        return f"TLTensorNetwork({', '.join(map(str, self.tensors))})"
    
    def __repr__(self):
        return str(self)
    
    @staticmethod
    def from_word(word: list[int], strands: int = None, normalize: bool = True) -> 'TLTensorNetwork':
        if strands is None:
            strands = max((abs(g) for g in word), default=0) + 1

        pos = list(range(strands))
        indices = [None] * strands
        n = strands
        tensors = []
        for g in word:
            idx = abs(g) - 1
            if g > 0:
                tensors.append(TLTensor.crossing([pos[idx], pos[idx + 1], n, n + 1], over=True))
            else:
                tensors.append(TLTensor.crossing([pos[idx], pos[idx + 1], n, n + 1], over=False))
            pos[idx] = n
            pos[idx+1] = n + 1
            indices[idx] = (len(tensors) - 1, 2)
            indices[idx+1] = (len(tensors) - 1, 3)
            n += 2

        for i in range(strands):
            if indices[i] is not None:
                idx, bidx = indices[i]
                bonds = tensors[idx].bonds
                bonds[bidx] = i
                tensors[idx] = tensors[idx].reindex(bonds)
            else:
                tensors.append(TLTensor([(TLDiagram.eye(1), [(1, 0)])], [i, i]))

        if normalize:
            writhe = sum(g // abs(g) for g in word)
            tensors.append(TLTensor.scalar([(1 if writhe % 2 == 0 else -1, -3 * writhe)]))

        return TLTensorNetwork(tensors)
    
    def _path_costs(self, path: list[tuple[int, int]]) -> list[int]:
        tensors = [(tensor.nterms, set(tensor.bonds)) for tensor in self.tensors]
        total = 0
        steps = []
        inter = 2
        for i, j in path:
            sa, ia = tensors.pop(max(i, j))
            sb, ib = tensors.pop(min(i, j))
            total += (len(ia) + len(ib)) * sa * sb
            steps.append(total)
            ic = ia.symmetric_difference(ib)
            sc = min(sa * sb, _catalan(len(ic) // 2))
            inter = max(inter, sc)
            tensors.append((sc, ic))
        return steps, inter

    def contract_info(self, optimize: str | cotengra.HyperOptimizer | list[tuple[int, int]] = "auto") -> ContractionInfo:
        if isinstance(optimize, str) and optimize == "sequential":
            # from_word stores crossings in braid order. A left fold keeps
            # a bounded frontier on narrow braids and needs no path search.
            # This remains an explicit choice: wide braids can fare poorly.
            count = len(self.tensors)
            path = [(0, 1)] if count > 1 else []
            path.extend((0, n - 1) for n in range(count - 1, 1, -1))
        else:
            tree = cotengra.array_contract_tree(
                inputs = [[str(idx) for idx in tensor.bonds] for tensor in self.tensors],
                shapes = [(2,)*len(tensor.bonds) for tensor in self.tensors],
                optimize=optimize
            )
            # Rust traces repeated bonds when merging tensors. Use the
            # tree binary path, without separate unary preprocessing steps.
            path = tree.get_path()
        steps, inter = self._path_costs(path)
        return ContractionInfo(math.log10(steps[-1]) if steps else 0.0, math.log2(inter), path)

    def contract(self, 
        poly: list[tuple[int, int]] = None, 
        optimize: str | cotengra.HyperOptimizer | list[tuple[int, int]] = "auto", 
        progbar: bool = False, 
        parallel: int = False
    ) -> 'TLTensorNetwork':
        info = self.contract_info(optimize)
        
        costs, _ = self._path_costs(info.path)
        prog = _tl_tensor.ProgressWrapper(costs[-1] if costs else 0, progbar)
        if parallel == "auto":
            pool = _tl_tensor.PoolWrapper(0)
        elif parallel == False:
            pool = _tl_tensor.PoolWrapper(-1)
        else:
            pool = _tl_tensor.PoolWrapper(parallel)
        if poly is None:
            poly = TLTensor._BRACKET_LOOP_POLY
        else:
            poly = _tl_tensor.PolyWrapper(*_poly_to_rust(poly))
        tensors = self.tensors.copy()
        for idx, (i, j) in enumerate(info.path):
            a = tensors.pop(max(i, j))
            b = tensors.pop(min(i, j))
            tensors.append(a._contract(b, prog, pool, poly))
            prog.set(costs[idx])
        prog.clear()
        return TLTensorNetwork(tensors)
    
    def contract_pair(self,
        i: int, j: int,
        poly: list[tuple[int, int]] = None, 
        progbar: bool = False,
        parallel: int = False
    ):
        a = self.tensors.pop(max(i, j))
        b = self.tensors.pop(min(i, j))
        self.tensors.append(a.contract(b, poly, progbar, parallel))

    def nterms(self) -> list[int]:
        return [tensor.nterms for tensor in self.tensors]
    
    def _proxy_qtn(self):
        import quimb.tensor as qtn
        import numpy as np
        return qtn.TensorNetwork([
            qtn.Tensor(np.empty((2,)*len(tensor.bonds)), inds=[str(i) for i in tensor.bonds])
            for tensor in self.tensors
        ])
    
    def draw(self):
        self._proxy_qtn().draw_3d_interactive()

    def _repr_html_(self):
        s = "<samp style='font-size: 12px;'>"
        s += "<details>"
        s += f"<summary>{_ach('TLTensorNetwork')}(tensors={_ach(len(self.tensors))})</summary>"
        for idx, tensor in enumerate(self.tensors):
            s += tensor._repr_html_()
            if idx >= 99:
                s += "<p>...</p>"
                break
        s += "</details></samp>"
        return s

def _catalan(n):
    return math.factorial(2*n) // (math.factorial(n + 1) * math.factorial(n))

def _tree_contraction_cost(tree: cotengra.ContractionTree) -> int:
    tensors = [(_catalan(len(bonds) // 2), set(bonds)) for bonds in tree.inputs]
    total = 0
    inter = 2
    for i, j in tree.get_path():
        sa, ia = tensors.pop(max(i, j))
        sb, ib = tensors.pop(min(i, j))
        total += (len(ia) + len(ib)) * sa * sb
        ic = ia.symmetric_difference(ib)
        sc = min(sa * sb, _catalan(len(ic) // 2))
        inter = max(inter, sc)
        tensors.append((sc, ic))
    return total, inter

class _TreeWrapper(cotengra.ContractionTree):
    _cost: float
    _mem: float
    def describe(self, info=None):
        return f"log10(cost) = {self._cost:.2f}, log2(mem) = {self._mem:.2f}"

def minimize_target(trial):
    stats = trial["tree"].contract_stats()
    if "flops" not in trial:
        trial["flops"] = stats["flops"]
    if "write" not in trial:
        trial["write"] = stats["write"]
    if "size" not in trial:
        trial["size"] = stats["size"]
    
    cost, mem = _tree_contraction_cost(trial['tree'])
    wrapped = _TreeWrapper.__new__(_TreeWrapper)
    wrapped.set_state_from(trial['tree'])
    wrapped._cost = math.log10(cost)
    wrapped._mem = math.log2(mem)
    trial['tree'] = wrapped
    return math.log10(cost)


class CotengraOptimizer(cotengra.HyperOptimizer):
    def __init__(self, parallel: int = False, progbar: bool = False, max_repeats: int = 128, **kwargs):
        kwargs['parallel'] = parallel
        kwargs['progbar'] = progbar
        kwargs['max_repeats'] = max_repeats
        kwargs['minimize'] = minimize_target
        # cotengra >= 0.8 refines every trial by subtree reconfiguration unless a refinement option
        # is given; that step needs a cotengra Objective (it optimises flops), not the TL cost above,
        # and fails every trial with it. None keeps the cotengra < 0.8 behaviour: no refinement.
        kwargs.setdefault('reconf_opts', None)

        if 'methods' not in kwargs:
            kwargs['methods'] = ['kahypar', 'greedy']
        super().__init__(**kwargs)

