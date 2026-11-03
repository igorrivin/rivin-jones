//! Direct TL braid updates on a frontier of at most twelve strands.
//!
//! A pairing of twenty-four endpoints fits in a u128 (five bits per partner).
//! For a bottom cup at i, (A I + A^-1 e_i) D = -A^-3 D; the inverse
//! crossing gives -A^3 D. Otherwise there are two shifted contributions.
//! Thus the total coefficient l1 norm grows by at most two per crossing.
//! Closing n strands multiplies it by at most 2^(n-1). The i128 path is
//! used only when c+n-1 <= 126, and every operation also checks overflow.
use crate::IntWrapper;
use pyo3::{exceptions::PyValueError, prelude::*};
use rug::Integer;
use std::collections::HashMap;

#[derive(Debug)]
struct Overflow;

trait Coefficient: Clone {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn add_scaled(&mut self, other: &Self, factor: i32) -> Result<(), Overflow>;
    fn to_python(&self) -> num_bigint::BigInt;
}

impl Coefficient for i128 {
    fn zero() -> Self { 0 }
    fn one() -> Self { 1 }
    fn is_zero(&self) -> bool { *self == 0 }
    fn add_scaled(&mut self, other: &Self, factor: i32) -> Result<(), Overflow> {
        let value = match factor {
            1 => self.checked_add(*other),
            -1 => self.checked_sub(*other),
            _ => other.checked_mul(factor as i128).and_then(|v| self.checked_add(v)),
        };
        *self = value.ok_or(Overflow)?;
        Ok(())
    }
    fn to_python(&self) -> num_bigint::BigInt { (*self).into() }
}

impl Coefficient for Integer {
    fn zero() -> Self { Integer::ZERO }
    fn one() -> Self { Integer::from(1) }
    fn is_zero(&self) -> bool { self == &0 }
    fn add_scaled(&mut self, other: &Self, factor: i32) -> Result<(), Overflow> {
        match factor {
            1 => *self += other,
            -1 => *self -= other,
            _ => *self += other * factor,
        }
        Ok(())
    }
    fn to_python(&self) -> num_bigint::BigInt {
        IntWrapper { inner: self.clone() }.to_numbigint()
    }
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 { (a, b) = (b, a % b); }
    a
}

#[derive(Clone)]
struct Polynomial<C> {
    min: isize,
    step: usize,
    coefficients: Vec<C>,
}

impl<C: Coefficient> Polynomial<C> {
    fn zero() -> Self { Self { min: 0, step: 0, coefficients: Vec::new() } }
    fn one() -> Self { Self { min: 0, step: 0, coefficients: vec![C::one()] } }
    fn is_zero(&self) -> bool { self.coefficients.is_empty() }
    fn max(&self) -> isize { self.min + ((self.coefficients.len()-1)*self.step) as isize }

    fn normalize(&mut self) {
        let Some(first) = self.coefficients.iter().position(|c| !c.is_zero()) else {
            *self = Self::zero();
            return;
        };
        let last = self.coefficients.iter().rposition(|c| !c.is_zero()).unwrap();
        self.coefficients.truncate(last+1);
        if first != 0 {
            self.min += (first*self.step) as isize;
            self.coefficients.drain(..first);
        }
        if self.coefficients.len() == 1 { self.step = 0; return; }
        let stride = self.coefficients.iter().enumerate()
            .filter(|(_, c)| !c.is_zero()).fold(0, |g, (i, _)| gcd(g, i));
        if stride > 1 {
            self.step *= stride;
            self.coefficients = std::mem::take(&mut self.coefficients)
                .into_iter().step_by(stride).collect();
        }
    }

    fn add_shifted(&mut self, other: &Self, shift: isize, factor: i32) -> Result<(), Overflow> {
        if other.is_zero() { return Ok(()); }
        if self.is_zero() && factor == 1 {
            *self = other.clone();
            self.min += shift;
            return Ok(());
        }
        let other_min = other.min + shift;
        let min = if self.is_zero() { other_min } else { self.min.min(other_min) };
        let max = if self.is_zero() { other.max()+shift } else { self.max().max(other.max()+shift) };
        let step = if self.is_zero() { other.step.max(1) }
                   else { gcd(gcd(self.step, other.step), self.min.abs_diff(other_min)).max(1) };
        let length = (max-min) as usize/step + 1;
        if self.is_zero() || self.step == 0 || self.step == step {
            let offset = if self.is_zero() { 0 } else { (self.min-min) as usize/step };
            self.coefficients.resize_with(length, C::zero);
            self.coefficients.rotate_right(offset);
        } else {
            let mut coefficients: Vec<C> = (0..length).map(|_| C::zero()).collect();
            for (i, c) in self.coefficients.drain(..).enumerate() {
                coefficients[(self.min-min) as usize/step+i*self.step/step] = c;
            }
            self.coefficients = coefficients;
        }
        let offset = (other_min-min) as usize/step;
        for (i, c) in other.coefficients.iter().enumerate() {
            self.coefficients[offset+i*other.step/step].add_scaled(c, factor)?;
        }
        self.min = min;
        self.step = step;
        self.normalize();
        Ok(())
    }

    fn to_python(&self) -> Vec<(num_bigint::BigInt, isize)> {
        self.coefficients.iter().enumerate().filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| (c.to_python(), self.min+(i*self.step) as isize)).collect()
    }
}

fn partner(diagram: u128, endpoint: usize) -> usize { ((diagram >> (5*endpoint)) & 31) as usize }
fn set_partner(diagram: &mut u128, endpoint: usize, other: usize) {
    *diagram = (*diagram & !(31_u128 << (5*endpoint))) | ((other as u128) << (5*endpoint));
}

struct Basis {
    diagrams: Vec<u128>,
    ids: HashMap<u128, usize>,
    actions: Vec<[usize; 11]>,
    strands: usize,
}

impl Basis {
    fn new(strands: usize) -> Self {
        let mut identity = 0;
        for i in 0..strands {
            set_partner(&mut identity, i, i+strands);
            set_partner(&mut identity, i+strands, i);
        }
        Self { diagrams: vec![identity], ids: HashMap::from([(identity, 0)]),
               actions: vec![[usize::MAX; 11]], strands }
    }

    fn action(&mut self, id: usize, generator: usize) -> usize {
        if self.actions[id][generator] != usize::MAX { return self.actions[id][generator]; }
        let mut diagram = self.diagrams[id];
        let left = self.strands+generator;
        let right = left+1;
        let a = partner(diagram, left);
        let b = partner(diagram, right);
        let target = if a == right { id } else {
            set_partner(&mut diagram, left, right);
            set_partner(&mut diagram, right, left);
            set_partner(&mut diagram, a, b);
            set_partner(&mut diagram, b, a);
            if let Some(&target) = self.ids.get(&diagram) { target } else {
                let target = self.diagrams.len();
                self.diagrams.push(diagram);
                self.ids.insert(diagram, target);
                self.actions.push([usize::MAX; 11]);
                target
            }
        };
        self.actions[id][generator] = target;
        target
    }

    fn closure_loops(&self, id: usize) -> usize {
        let mut seen = 0_u32;
        let mut loops = 0;
        for start in 0..2*self.strands {
            if seen & (1 << start) != 0 { continue; }
            loops += 1;
            let mut current = start;
            while seen & (1 << current) == 0 {
                let next = partner(self.diagrams[id], current);
                seen |= (1 << current) | (1 << next);
                current = if next < self.strands { next+self.strands } else { next-self.strands };
            }
        }
        loops
    }
}

fn compute<C: Coefficient>(word: &[isize], strands: usize) -> Result<Polynomial<C>, Overflow> {
    let mut basis = Basis::new(strands);
    let mut current = vec![Polynomial::<C>::one()];
    for &g in word {
        let generator = g.unsigned_abs()-1;
        let sign = g.signum();
        let mut next: Vec<Polynomial<C>> = (0..basis.diagrams.len()).map(|_| Polynomial::zero()).collect();
        for (id, mut polynomial) in current.into_iter().enumerate() {
            if polynomial.is_zero() { continue; }
            let target = basis.action(id, generator);
            next.resize_with(basis.diagrams.len(), Polynomial::zero);
            if target == id {
                next[id].add_shifted(&polynomial, -3*sign, -1)?;
            } else {
                next[target].add_shifted(&polynomial, -sign, 1)?;
                polynomial.min += sign;
                // All e_i outputs have a cup at i; this identity output does
                // not, so it cannot collide with another contribution.
                debug_assert!(next[id].is_zero());
                next[id] = polynomial;
            }
        }
        current = next;
    }
    let writhe: isize = word.iter().map(|g| g.signum()).sum();
    let mut result = Polynomial::<C>::zero();
    for (id, p) in current.iter().enumerate() {
        if p.is_zero() { continue; }
        let power = basis.closure_loops(id)-1;
        // delta^power = (-1)^power sum_j binom(power,j) A^(2power-4j).
        let sign = if (writhe + power as isize) % 2 == 0 { 1 } else { -1 };
        let mut binomial = 1;
        for j in 0..=power {
            result.add_shifted(p, 2*power as isize-4*j as isize-3*writhe, sign*binomial)?;
            if j < power { binomial = binomial*(power-j) as i32/(j+1) as i32; }
        }
    }
    Ok(result)
}

#[pyfunction]
#[pyo3(signature=(word, strands, arithmetic="auto"))]
pub(crate) fn braid_jones(word: Vec<isize>, strands: usize, arithmetic: &str)
    -> PyResult<(Vec<(num_bigint::BigInt, isize)>, &'static str)> {
    if !(1..=12).contains(&strands) {
        return Err(PyValueError::new_err("frontier computation requires 1 through 12 strands"));
    }
    if word.iter().any(|g| *g == 0 || g.unsigned_abs() >= strands) {
        return Err(PyValueError::new_err("braid generators must satisfy 1 <= abs(g) < strands"));
    }
    if arithmetic != "auto" && arithmetic != "gmp" {
        return Err(PyValueError::new_err("arithmetic must be auto or gmp"));
    }
    if arithmetic == "auto" && strands == 2 {
        // B_2 is infinite cyclic. Sum the signed generators and evaluate the
        // exact closed form, linear in input plus output size. For m>0,
        // J = (-1)^(m+1) A^(-2m+2)
        //     + sum_{k=1}^{m-1} (-1)^(m+k) A^(-2m-2-4k).
        // This follows by summing (A I + A^-1 e)^m with e^2=delta e.
        // A negative exponent mirrors A; exponent zero is the two-component unlink.
        let writhe: isize = word.iter().sum();
        if writhe == 0 {
            return Ok((vec![(num_bigint::BigInt::from(-1), -2),
                            (num_bigint::BigInt::from(-1), 2)], "b2-formula"));
        }
        let m = writhe.abs();
        let sign = writhe.signum();
        let mut coefficients = Vec::with_capacity(m as usize);
        coefficients.push((num_bigint::BigInt::from(if m % 2 == 0 { -1 } else { 1 }), sign*(-2*m+2)));
        for k in 1..m {
            coefficients.push((num_bigint::BigInt::from(if (m+k) % 2 == 0 { 1 } else { -1 }),
                               sign*(-2*m-2-4*k)));
        }
        return Ok((coefficients, "b2-formula"));
    }
    if arithmetic == "auto" && word.len().checked_add(strands-1).is_some_and(|b| b <= 126) {
        if let Ok(result) = compute::<i128>(&word, strands) {
            return Ok((result.to_python(), "i128"));
        }
        // Defensive fallback: an overflow can never return a wrapped answer.
    }
    let result = compute::<Integer>(&word, strands).expect("GMP integer addition cannot overflow");
    Ok((result.to_python(), "gmp"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_integer_overflow() {
        let mut maximum = i128::MAX;
        assert!(maximum.add_scaled(&1, 1).is_err());
        let mut minimum = i128::MIN;
        assert!(minimum.add_scaled(&1, -1).is_err());
        assert!(0_i128.add_scaled(&i128::MAX, 2).is_err());
    }

    #[test]
    fn integer_modes_and_signed_lattice_addition() {
        for trial in 0..300 {
            let n = 2+trial % 6;
            let word: Vec<isize> = (0..(trial % 25)).map(|i| {
                let sign = if (i+trial) % 3 == 0 { -1 } else { 1 };
                sign*(1+(i*7+trial) % (n-1)) as isize
            }).collect();
            assert_eq!(compute::<i128>(&word, n).unwrap().to_python(),
                       compute::<Integer>(&word, n).unwrap().to_python());
        }
        let mut small = Polynomial::<i128>::zero();
        let mut big = Polynomial::<Integer>::zero();
        for i in 0..120 {
            let factor = if i % 3 == 0 { -7 } else { 5 };
            small.add_shifted(&Polynomial::one(), (i % 17)-8, factor).unwrap();
            big.add_shifted(&Polynomial::one(), (i % 17)-8, factor).unwrap();
            assert_eq!(small.to_python(), big.to_python());
        }
    }

    #[test]
    fn shifted_sums_against_independent_integer_map() {
        use std::collections::BTreeMap;
        let mut state = 771_u64;
        for trial in 0..200 {
            let mut result = Polynomial::<i128>::zero();
            let mut expected = BTreeMap::<isize, num_bigint::BigInt>::new();
            for term in 0..12 {
                let min = (term*13+trial) as isize % 31-15;
                let step = 1+(term+trial) % 7;
                let shift = (term*3) as isize-20;
                let factor = if term % 2 == 0 { -7 } else { 5 };
                let mut p = Polynomial { min, step, coefficients: Vec::new() };
                for i in 0..1+(trial+term) % 23 {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let c = (state as i64 as i128) << (trial % 30);
                    p.coefficients.push(c);
                    *expected.entry(min+shift+(i*step) as isize).or_default()
                        += num_bigint::BigInt::from(c)*factor;
                }
                p.normalize();
                result.add_shifted(&p, shift, factor).unwrap();
            }
            let actual = result.to_python().into_iter().map(|(c, e)| (e, c)).collect::<BTreeMap<_, _>>();
            expected.retain(|_, c| c != &num_bigint::BigInt::from(0));
            assert_eq!(actual, expected);
        }
    }
}
