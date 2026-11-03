//! Exact integer Laurent polynomials on their actual exponent lattice.
//!
//! min + i*step indexes the coefficients. The stride is derived from support,
//! not assumed from a Jones convention. Addition refines the lattice when
//! residues differ. Large products use signed, bounded Kronecker substitution.
use crate::IntWrapper;
use num_traits::Zero;
use rug::{integer::Order, ops::NegAssign, Integer};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Laurent {
    min: isize,
    step: usize,
    coefficients: Vec<IntWrapper>,
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 { (a, b) = (b, a % b); }
    a
}

impl Laurent {
    pub(crate) fn new(min: isize, coefficients: Vec<IntWrapper>) -> Self {
        Self::from_parts(min, 1, coefficients)
    }

    fn from_parts(min: isize, step: usize, coefficients: Vec<IntWrapper>) -> Self {
        let mut p = Self { min, step, coefficients };
        p.normalize();
        p
    }

    fn normalize(&mut self) {
        let Some(first) = self.coefficients.iter().position(|c| !c.is_zero()) else {
            self.min = 0;
            self.step = 0;
            self.coefficients.clear();
            return;
        };
        let last = self.coefficients.iter().rposition(|c| !c.is_zero()).unwrap();
        self.coefficients.truncate(last + 1);
        if first != 0 {
            self.min += (first * self.step) as isize;
            self.coefficients.drain(..first);
        }
        if self.coefficients.len() == 1 {
            self.step = 0;
            return;
        }
        let stride = self.coefficients.iter().enumerate()
            .filter(|(_, c)| !c.is_zero()).fold(0, |g, (i, _)| gcd(g, i));
        if stride > 1 {
            self.step *= stride;
            self.coefficients = std::mem::take(&mut self.coefficients)
                .into_iter().step_by(stride).collect();
        }
    }

    pub(crate) fn is_empty(&self) -> bool { self.coefficients.is_empty() }

    pub(crate) fn iter(&self) -> impl Iterator<Item=(isize, &IntWrapper)> {
        self.coefficients.iter().enumerate()
            .map(move |(i, c)| (self.min + (i*self.step) as isize, c))
    }

    fn max(&self) -> isize {
        self.min + ((self.coefficients.len()-1)*self.step) as isize
    }
}

impl std::ops::AddAssign<&Laurent> for Laurent {
    fn add_assign(&mut self, other: &Laurent) {
        if other.is_empty() { return; }
        if self.is_empty() { *self = other.clone(); return; }
        let min = self.min.min(other.min);
        let max = self.max().max(other.max());
        let step = gcd(gcd(self.step, other.step), self.min.abs_diff(other.min)).max(1);
        let length = (max-min) as usize / step + 1;
        let left_offset = (self.min-min) as usize / step;
        let right_offset = (other.min-min) as usize / step;
        if (self.step == 0 || self.step == step) && (other.step == 0 || other.step == step) {
            self.coefficients.resize(length, IntWrapper::zero());
            // The new trailing slots include at least left_offset zeros.
            self.coefficients.rotate_right(left_offset);
            for (i, c) in other.coefficients.iter().enumerate() {
                self.coefficients[right_offset+i].inner += &c.inner;
            }
        } else {
            let mut coefficients = vec![IntWrapper::zero(); length];
            for (e, c) in self.iter().chain(other.iter()) {
                coefficients[(e-min) as usize / step].inner += &c.inner;
            }
            self.coefficients = coefficients;
        }
        self.min = min;
        self.step = step;
        self.normalize();
    }
}

fn scale_shift(p: &Laurent, scalar: &IntWrapper, shift: isize) -> Laurent {
    let mut result = p.clone();
    result.min += shift;
    if scalar.inner == -1 {
        for c in &mut result.coefficients { c.inner.neg_assign(); }
    } else if scalar.inner != 1 {
        for c in &mut result.coefficients { c.inner *= &scalar.inner; }
    }
    result
}

pub(crate) fn multiply(a: &Laurent, b: &Laurent) -> Laurent {
    if a.is_empty() || b.is_empty() { return Laurent::new(0, Vec::new()); }
    if a.coefficients.len() == 1 { return scale_shift(b, &a.coefficients[0], a.min); }
    if b.coefficients.len() == 1 { return scale_shift(a, &b.coefficients[0], b.min); }
    let av: Vec<_> = a.iter().filter(|(_, c)| !c.is_zero()).collect();
    let bv: Vec<_> = b.iter().filter(|(_, c)| !c.is_zero()).collect();
    let step = gcd(a.step, b.step).max(1);
    let length = ((a.max()-a.min) + (b.max()-b.min)) as usize / step + 1;
    if av.len().min(bv.len()) <= 4 || av.len().saturating_mul(bv.len()) <= 256 {
        let mut result = vec![IntWrapper::zero(); length];
        for &(ae, ac) in &av {
            for &(be, bc) in &bv {
                let index = ((ae-a.min) + (be-b.min)) as usize / step;
                result[index].inner += &ac.inner * &bc.inner;
            }
        }
        return Laurent::from_parts(a.min+b.min, step, result);
    }
    kronecker(a, b, &av, &bv, step)
}

/// Each coefficient has magnitude < r*2^(u+v), where r is the smaller
/// nonzero count and u,v are the largest input coefficient bit lengths.
/// Byte-align k >= u+v+ceil(log2(r))+1. The product coefficients are strictly
/// between -B/2 and B/2 for B=2^k, so balanced-base decoding is unique.
fn kronecker(a: &Laurent, b: &Laurent, av: &[(isize, &IntWrapper)],
             bv: &[(isize, &IntWrapper)], step: usize) -> Laurent {
    let abits = av.iter().map(|(_, c)| c.inner.significant_bits()).max().unwrap() as usize;
    let bbits = bv.iter().map(|(_, c)| c.inner.significant_bits()).max().unwrap() as usize;
    let r = av.len().min(bv.len());
    let overlap_bits = (usize::BITS - (r-1).leading_zeros()) as usize;
    let width = (abits + bbits + overlap_bits + 1).div_ceil(8);
    let bits = width * 8;

    fn pack(values: &[(isize, &IntWrapper)], min: isize, count: usize,
            step: usize, width: usize) -> Integer {
        let mut positive = vec![0_u8; count*width];
        let mut negative = vec![0_u8; count*width];
        for &(e, c) in values {
            let offset = (e-min) as usize / step * width;
            let buffer = if c.inner.is_negative() { &mut negative } else { &mut positive };
            c.inner.write_digits(&mut buffer[offset..offset+width], Order::Lsf);
        }
        Integer::from_digits(&positive, Order::Lsf) - Integer::from_digits(&negative, Order::Lsf)
    }

    let alen = (a.max()-a.min) as usize / step + 1;
    let blen = (b.max()-b.min) as usize / step + 1;
    let value = pack(av, a.min, alen, step, width) * pack(bv, b.min, blen, step, width);
    let negative = value.is_negative();
    let digits = value.to_digits::<u8>(Order::Lsf);
    let base = Integer::from(1) << bits;
    let half = Integer::from(1) << (bits-1);
    let mut result = Vec::with_capacity(alen+blen-1);
    let mut carry = false;
    for i in 0..alen+blen-1 {
        let start = i*width;
        let end = (start+width).min(digits.len());
        let mut c = if start < end {
            Integer::from_digits(&digits[start..end], Order::Lsf)
        } else { Integer::ZERO };
        if carry { c += 1; }
        carry = c >= half;
        if carry { c -= &base; }
        if negative { c.neg_assign(); }
        result.push(IntWrapper { inner: c });
    }
    debug_assert!(!carry, "Kronecker coefficient bound was insufficient");
    Laurent::from_parts(a.min+b.min, step, result)
}

#[cfg(test)]
mod tests {
    use super::*;
    type Reference = series::Polynomial<IntWrapper>;

    fn dense(p: &Laurent) -> Reference {
        if p.is_empty() { return Reference::new(0, Vec::new()); }
        let mut coefficients = vec![IntWrapper::zero(); (p.max()-p.min) as usize+1];
        for (e, c) in p.iter() { coefficients[(e-p.min) as usize] = c.clone(); }
        Reference::new(p.min, coefficients)
    }

    fn check(actual: &Laurent, expected: &Reference) {
        let a: Vec<_> = actual.iter().filter(|(_, c)| !c.is_zero()).collect();
        let b: Vec<_> = expected.iter().filter(|(_, c)| !c.is_zero()).collect();
        assert_eq!(a, b);
    }

    #[test]
    fn zero_monomials_and_cancellation() {
        let zero = Laurent::new(0, Vec::new());
        let a = Laurent::new(-9, vec![1,0,0,0,-3,0,0,0,2].into_iter().map(IntWrapper::from).collect());
        assert_eq!(a.step, 4);
        assert_eq!(multiply(&zero, &a), zero);
        for c in [-7, -1, 1, 5] {
            let b = Laurent::new(12, vec![IntWrapper::from(c)]);
            check(&multiply(&a, &b), &(&dense(&a)*&dense(&b)));
        }
        let mut cancelled = a.clone();
        cancelled += &scale_shift(&a, &IntWrapper::from(-1), 0);
        assert_eq!(cancelled, zero);
    }

    #[test]
    fn support_residues_signed_products_and_sums() {
        let mut state = 771_u64;
        for trial in 0..500 {
            let mut make = |length: usize, step: usize, min: isize| {
                let mut coefficients = vec![IntWrapper::zero(); length*step];
                for i in 0..length {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let mut coefficient = Integer::from((state % 101) as i32 - 50);
                    coefficient <<= trial % 257;
                    coefficients[i*step] = IntWrapper { inner: coefficient };
                }
                Laurent::new(min, coefficients)
            };
            let a = make(1+trial % 43, 1+trial % 7, -21);
            let b = make(1+trial % 37, 1+(trial/7) % 7, trial as isize % 31-17);
            let expected = &dense(&a)*&dense(&b);
            check(&multiply(&a, &b), &expected);
            let mut sum = a.clone();
            sum += &b;
            check(&sum, &(&dense(&a)+&dense(&b)));
            if !a.is_empty() && !b.is_empty() {
                let av: Vec<_> = a.iter().filter(|(_, c)| !c.is_zero()).collect();
                let bv: Vec<_> = b.iter().filter(|(_, c)| !c.is_zero()).collect();
                check(&kronecker(&a, &b, &av, &bv, gcd(a.step, b.step).max(1)), &expected);
            }
        }
    }
}
