use std::collections::HashMap;
use pyo3::prelude::*;
use num_traits::Zero;

type LPoly = laurent::Laurent;
mod laurent;
mod frontier;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct IntWrapper {
    inner: rug::Integer
}

impl IntWrapper {
    fn from_numbigint(val: &num_bigint::BigInt) -> Self {
        use rug::ops::NegAssign;

        let (sign, digits) = val.to_bytes_be();
        let mut inner = rug::Integer::from_digits(&digits, rug::integer::Order::Msf);
        if let num_bigint::Sign::Minus = sign {
            inner.neg_assign();
        }

        IntWrapper {
            inner
        }
    }

    fn to_numbigint(&self) -> num_bigint::BigInt {
        let bytes = self.inner.to_digits::<u8>(rug::integer::Order::Msf);
        let sign = if self.inner.is_positive() {
            num_bigint::Sign::Plus
        } else if self.inner.is_zero() {
            num_bigint::Sign::NoSign
        } else {
            num_bigint::Sign::Minus
        };
        num_bigint::BigInt::from_bytes_be(sign, &bytes)
    }
}

impl From<isize> for IntWrapper {
    fn from(value: isize) -> Self {
        IntWrapper { inner: rug::Integer::from(value) }
    }
}

impl std::ops::Add<IntWrapper> for IntWrapper {
    type Output = IntWrapper;
    fn add(self, other: Self) -> Self {
        IntWrapper { inner: self.inner + other.inner }
    }
}

impl std::ops::AddAssign<&IntWrapper> for IntWrapper {
    fn add_assign(&mut self, rhs: &IntWrapper) {
        self.inner += &rhs.inner;
    }
}

impl std::ops::AddAssign<IntWrapper> for IntWrapper {
    fn add_assign(&mut self, rhs: IntWrapper) {
        self.inner += rhs.inner;
    }
}

impl std::ops::SubAssign<&IntWrapper> for IntWrapper {
    fn sub_assign(&mut self, rhs: &IntWrapper) {
        self.inner -= &rhs.inner;
    }
}

impl std::ops::SubAssign<IntWrapper> for IntWrapper {
    fn sub_assign(&mut self, rhs: IntWrapper) {
        self.inner -= rhs.inner;
    }
}

impl num_traits::Zero for IntWrapper {
    fn zero() -> Self {
        IntWrapper { inner: rug::Integer::ZERO }
    }

    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }
}

impl std::ops::Mul<IntWrapper> for IntWrapper {
    type Output = IntWrapper;
    fn mul(self, rhs: IntWrapper) -> Self::Output {
        IntWrapper { inner: self.inner * rhs.inner }
    }
}

impl std::ops::Mul<IntWrapper> for &IntWrapper {
    type Output = IntWrapper;
    fn mul(self, rhs: IntWrapper) -> Self::Output {
        IntWrapper { inner: &self.inner * rhs.inner }
    }
}

impl std::ops::Mul<&IntWrapper> for &IntWrapper {
    type Output = IntWrapper;
    fn mul(self, rhs: &IntWrapper) -> Self::Output {
        IntWrapper { inner: (&self.inner * &rhs.inner).into() }
    }
}

impl std::ops::MulAssign<&IntWrapper> for IntWrapper {
    fn mul_assign(&mut self, rhs: &IntWrapper) {
        self.inner *= &rhs.inner;
    }
}

impl std::ops::MulAssign<IntWrapper> for IntWrapper {
    fn mul_assign(&mut self, rhs: IntWrapper) {
        self.inner *= rhs.inner;
    }
}

impl num_traits::One for IntWrapper {
    fn one() -> Self {
        IntWrapper { inner: rug::Integer::ONE.clone() }
    }

    fn is_one(&self) -> bool {
        self.inner.is_one()
    }
}

impl std::ops::Neg for IntWrapper {
    type Output = IntWrapper;
    fn neg(self) -> Self::Output {
        IntWrapper { inner: -self.inner }
    }
}

impl std::ops::Neg for &IntWrapper {
    type Output = IntWrapper;
    fn neg(self) -> Self::Output {
        IntWrapper { inner: (-&self.inner).into() }
    }
}

#[pyclass]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct TLDiagramInner {
    #[pyo3(get, set)]
    idx: Vec<u8>,
    #[pyo3(get, set)]
    loops: usize
}

#[pymethods]
impl TLDiagramInner {
    #[new]
    fn new(idx: Vec<u8>, loops: usize) -> TLDiagramInner {
        TLDiagramInner { idx, loops }
    }
    
    fn tensor(&self, other: &TLDiagramInner) -> TLDiagramInner {
        let sn = self.idx.len() / 2;
        let on = other.idx.len() / 2;
        let mut new_idx = vec![0; self.idx.len() + other.idx.len()];

        for mut i in 0..self.idx.len() {
            let mut idx = self.idx[i] as usize;
            if i >= sn { i += on; }
            if idx >= sn { idx += on; }
            new_idx[i] = idx as u8;
            new_idx[idx] = i as u8;
        }

        for mut i in 0..other.idx.len() {
            let mut idx = other.idx[i] as usize;
            if i >= on { i += sn; }
            i += sn;
            if idx >= on { idx += sn; }
            idx += sn;
            new_idx[i] = idx as u8;
            new_idx[idx] = i as u8;
        }

        TLDiagramInner {
            idx: new_idx, loops: self.loops + other.loops
        }
    }

    fn trace(&mut self, u: usize, v: usize) {
        let a = self.idx[u] as usize;
        let b = self.idx[v] as usize;
        if a == v && b == u {
            self.loops += 1;
        }
        self.idx[a] = b as u8;
        self.idx[b] = a as u8;

        for i in 0..self.idx.len() {
            let val = self.idx[i] as usize;
            if val > u { self.idx[i] -= 1; }
            if val > v { self.idx[i] -= 1; }
        }

        self.idx.remove(u.max(v));
        self.idx.remove(u.min(v));
    }
}

#[pyclass]
struct PoolWrapper {
    pool: Option<rayon::ThreadPool>
}

#[pymethods]
impl PoolWrapper {
    #[new]
    fn new(threads: isize) -> PoolWrapper {
        if threads < 0 {
            PoolWrapper { pool: None }
        } else if threads == 0 {
            PoolWrapper { pool: Some(rayon::ThreadPoolBuilder::new().num_threads(num_cpus::get_physical()).build().unwrap()) }
        } else {
            PoolWrapper { pool: Some(rayon::ThreadPoolBuilder::new().num_threads(threads as usize).build().unwrap()) }
        }
    }
}

impl PoolWrapper {
    fn par_for_each<I: IntoIterator + Send>(&self, iter: I, f: impl Fn(<I as IntoIterator>::Item) + Send + Sync)
        where I: rayon::iter::IntoParallelIterator<Item=<I as IntoIterator>::Item> {
        use rayon::iter::ParallelIterator;

        if let Some(pool) = &self.pool {
            pool.install(|| {
                iter.into_par_iter()
                    .for_each(f)
            })
        } else {
            iter.into_iter()
                .for_each(f)
        }
    }
}

#[pyclass]
struct ProgressWrapper {
    bar: Option<indicatif::ProgressBar>
}

#[pymethods]
impl ProgressWrapper {
    #[new]
    fn new(total: usize, show: bool) -> ProgressWrapper {
        ProgressWrapper {
            bar: if show {
                let prog = indicatif::ProgressBar::new(total as u64);
                prog.set_style(indicatif::ProgressStyle::with_template("[{elapsed}] {wide_bar:.cyan/blue} {percent}% {pos:>7}/{len:7} [{eta} eta, {decimal_bytes_per_sec}]")
                    .unwrap()
                    .progress_chars("##-"));
                Some(prog)
            } else {
                None
            }
        }
    }

    fn inc(&self) {
        if let Some(bar) = &self.bar {
            bar.inc(1);
        }
    }

    fn set(&self, value: usize) {
        if let Some(bar) = &self.bar {
            bar.set_position(value as u64)
        }
    }

    fn clear(&self) {
        if let Some(bar) = &self.bar {
            bar.finish();
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct PolyWrapper {
    inner: LPoly
}

#[pymethods]
impl PolyWrapper {
    #[new]
    fn new(min: isize, coeffs: Vec<num_bigint::BigInt>) -> PyResult<Self> {
        if coeffs.len() > 1000001 || min.checked_add(coeffs.len().saturating_sub(1) as isize).is_none() {
            return Err(pyo3::exceptions::PyValueError::new_err("polynomial input span or exponent exceeds the native arithmetic guard"));
        }
        Ok(PolyWrapper {
            inner: LPoly::new(min, coeffs.iter().map(IntWrapper::from_numbigint).collect())
        })
    }

    fn __add__(&self, other: &Self) -> PyResult<Self> {
        check_polynomial_span(&self.inner, &other.inner, false)?;
        let mut inner = self.inner.clone();
        inner += &other.inner;
        Ok(Self { inner })
    }

    fn __sub__(&self, other: &Self) -> PyResult<Self> {
        self.__add__(&other.__neg__())
    }

    fn __mul__(&self, other: &Self) -> PyResult<Self> {
        check_polynomial_span(&self.inner, &other.inner, true)?;
        Ok(Self { inner: laurent::multiply(&self.inner, &other.inner) })
    }

    fn __neg__(&self) -> Self {
        Self { inner: laurent::multiply(&self.inner, &LPoly::new(0, vec![IntWrapper::from(-1)])) }
    }

    fn __eq__(&self, other: &Self) -> bool { self.inner == other.inner }
    fn __bool__(&self) -> bool { !self.inner.is_empty() }

    fn __pow__(&self, exponent: isize, modulo: Option<isize>) -> PyResult<Self> {
        use pyo3::exceptions::PyValueError;
        if modulo.is_some() { return Err(PyValueError::new_err("modular polynomial power is unsupported")); }
        let terms: Vec<_> = self.inner.iter().filter(|(_, c)| !c.is_zero()).collect();
        if terms.len() == 1 && (terms[0].1.inner == 1 || terms[0].1.inner == -1) {
            let power = terms[0].0.checked_mul(exponent)
                .ok_or_else(|| PyValueError::new_err("polynomial exponent overflow"))?;
            let sign = if terms[0].1.inner == -1 && exponent % 2 != 0 { -1 } else { 1 };
            return Ok(Self { inner: LPoly::new(power, vec![IntWrapper::from(sign)]) });
        }
        if !(0..=100000).contains(&exponent) {
            return Err(PyValueError::new_err("power requires a nonnegative bounded exponent or a unit monomial"));
        }
        let mut result = Self { inner: LPoly::new(0, vec![IntWrapper::from(1)]) };
        let mut base = self.clone();
        let mut n = exponent;
        while n != 0 {
            if n % 2 != 0 { result = result.__mul__(&base)?; }
            n /= 2;
            if n != 0 { base = base.__mul__(&base)?; }
        }
        Ok(result)
    }

    #[getter]
    fn terms(&self) -> Vec<(num_bigint::BigInt, isize)> {
        self.inner.iter().filter(|(_, c)| !c.is_zero())
            .map(|(e, c)| (c.to_numbigint(), e)).collect()
    }
}

fn check_polynomial_span(a: &LPoly, b: &LPoly, multiply: bool) -> PyResult<()> {
    use pyo3::exceptions::PyValueError;
    if a.is_empty() || b.is_empty() { return Ok(()); }
    let mut ai = a.iter();
    let mut bi = b.iter();
    let amin = ai.next().unwrap().0;
    let bmin = bi.next().unwrap().0;
    let amax = ai.last().map_or(amin, |(e, _)| e);
    let bmax = bi.last().map_or(bmin, |(e, _)| e);
    let (low, high) = if multiply { (amin.checked_add(bmin), amax.checked_add(bmax)) }
                      else { (Some(amin.min(bmin)), Some(amax.max(bmax))) };
    let span = low.zip(high).and_then(|(l, h)| h.checked_sub(l));
    if span.is_none_or(|s| s > 1000000) {
        return Err(PyValueError::new_err("exact polynomial span exceeds the native arithmetic guard"));
    }
    Ok(())
}

#[pyclass]
#[derive(Debug, Clone)]
struct TLTensorInner {
    terms: Vec<(TLDiagramInner, LPoly)>,
    #[pyo3(get, set)]
    bonds: Vec<usize>
}

#[pymethods]
impl TLTensorInner {
    #[new]
    fn new(terms: Vec<(TLDiagramInner, isize, Vec<num_bigint::BigInt>)>, bonds: Vec<usize>) -> TLTensorInner {
        TLTensorInner {
            bonds, terms: terms.into_iter()
                .map(|(d, m, c)| (d, m, c.iter()
                    .map(|i| IntWrapper::from_numbigint(&i))
                    .collect()
                ))
                .map(|(d, m, c)| (d, LPoly::new(m, c)))
                .filter(|(_, p)| !p.is_empty())
                .collect()
        }
    }

    #[getter]
    fn terms(&self) -> Vec<(TLDiagramInner, Vec<(num_bigint::BigInt, isize)>)> {
        self.terms.iter()
            .map(|(d, p)| (d.clone(), p.iter()
                .map(|(a, b)| (b.to_numbigint(), a))
                .filter(|(a, _)| !a.is_zero())
                .collect::<Vec<_>>()))
            .collect()
    }

    #[getter]
    fn nterms(&self) -> usize {
        return self.terms.len()
    }

    fn contract_plan(&self, other: &TLTensorInner) -> (Vec<(usize, usize)>, Vec<usize>) {
        let sn = self.bonds.len() / 2;
        let on = other.bonds.len() / 2;
        let mut cbonds = Vec::new();
        cbonds.extend_from_slice(&self.bonds[..sn]);
        cbonds.extend_from_slice(&other.bonds[..on]);
        cbonds.extend_from_slice(&self.bonds[sn..]);
        cbonds.extend_from_slice(&other.bonds[on..]);

        let mut plan = Vec::new();
        'outer: loop {
            if cbonds.is_empty() {
                break
            }

            for (i, &b) in cbonds[..cbonds.len() - 1].iter().enumerate() {
                let Some(jp) = cbonds[i+1..].iter()
                    .position(|&bp| bp == b) else { continue };
                let j = i + 1 + jp;

                if j < cbonds.len() - 1 && cbonds[j+1..].iter().any(|&bp| bp == b) {
                    panic!("Hyperbonds not supported");
                }

                if (i < cbonds.len() / 2 && j < cbonds.len() / 2) 
                    || (i >= cbonds.len() / 2 && j >= cbonds.len() / 2) {
                    panic!("Cannot trace indices that are on the same side of a diagram");
                } 

                plan.push((i, j));
                cbonds.remove(j);
                cbonds.remove(i);
                continue 'outer;
            }
            break
        }

        (plan, cbonds)
    }

    fn add(&self, other: &TLTensorInner) -> TLTensorInner {
        assert_eq!(self.bonds, other.bonds);

        let mut cterms = HashMap::new();
        for (d, p) in &self.terms {
            cterms.entry(d.clone())
                .and_modify(|pc| *pc += p)
                .or_insert(p.clone());
        }
        for (d, p) in &other.terms {
            cterms.entry(d.clone())
                .and_modify(|pc| *pc += p)
                .or_insert(p.clone());
        }
        TLTensorInner { terms: cterms.into_iter().collect(), bonds: self.bonds.clone() }
    }

    fn contract(&self, other: &TLTensorInner, prog: &mut ProgressWrapper, pool: &mut PoolWrapper, loop_poly: &PolyWrapper) -> TLTensorInner {
        let (plan, cbonds) = self.contract_plan(other);

        let product = |idx: usize| {
            let i = idx / other.terms.len();
            let j = idx % other.terms.len();
            let (da, pa) = &self.terms[i];
            let (db, pb) = &other.terms[j];
            let mut dc = da.tensor(db);
            let mut pc = laurent::multiply(pa, pb);
            for &(i, j) in &plan {
                dc.trace(i, j);
            }

            let nloops = if dc.idx.len() > 0 {
                dc.loops
            } else {
                assert!(dc.loops >= 1);
                dc.loops - 1
            };

            for _ in 0..nloops {
                dc.loops -= 1;
                pc = laurent::multiply(&pc, &loop_poly.inner);
            }

            prog.inc();
            (dc, pc)
        };

        let count = self.terms.len()*other.terms.len();
        let terms = if pool.pool.is_none() {
            // Serial contraction needs no shard locks or concurrent-map setup.
            let mut cterms = HashMap::new();
            for idx in 0..count {
                let (dc, pc) = product(idx);
                cterms.entry(dc).and_modify(|p| *p += &pc).or_insert(pc);
            }
            cterms.into_iter().filter(|(_, p)| !p.is_empty()).collect()
        } else {
            let cterms = dashmap::DashMap::new();
            pool.par_for_each(0..count, |idx| {
                let (dc, pc) = product(idx);
                cterms.entry(dc).and_modify(|p| *p += &pc).or_insert(pc);
            });
            cterms.into_iter().filter(|(_, p)| !p.is_empty()).collect()
        };
        TLTensorInner { terms, bonds: cbonds }
    }
}

#[pymodule]
fn _tl_tensor(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(frontier::braid_jones, m)?)?;
    m.add_class::<TLDiagramInner>()?;
    m.add_class::<TLTensorInner>()?;
    m.add_class::<ProgressWrapper>()?;
    m.add_class::<PoolWrapper>()?;
    m.add_class::<PolyWrapper>()?;
    Ok(())
}
