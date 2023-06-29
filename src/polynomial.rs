use std::{
    iter::Sum,
    ops::{Add, Mul, MulAssign},
};

use bls12_381::Scalar;

#[derive(Copy, Clone)]
pub struct Polynomial<const TERMS: usize> {
    terms: [Scalar; TERMS],
}

impl<const TERMS: usize> Polynomial<TERMS> {
    /// The polynomial t(x) = (x - w_1)(x - w_2)...
    /// for all roots of unity. The polynomial with 2^k + 1 terms
    /// will have its highest degree at 2^k as expected.
    pub fn t_x(degree: usize) -> Self {
        let mut poly = Self::zero();
        poly.terms[0] = -Scalar::one();
        poly.terms[degree] = Scalar::one();
        poly
    }

    pub const fn zero() -> Self {
        Polynomial {
            terms: [Scalar::zero(); TERMS],
        }
    }

    pub fn identity() -> Self {
        let mut poly = Polynomial {
            terms: [Scalar::zero(); TERMS],
        };
        poly.terms[0] = Scalar::one();
        poly
    }

    /// Get the polynomial corresponding to (x - `const_term`).
    pub fn x_minus_const(const_term: Scalar) -> Self {
        let mut poly = Polynomial {
            terms: [Scalar::zero(); TERMS],
        };
        poly.terms[0] = -const_term;
        poly.terms[1] = Scalar::one();
        poly
    }

    /// Get the degree of this polynomial. Returns usize::MAX for zero-polynomials
    pub fn degree(&self) -> usize {
        self.terms
            .iter()
            .enumerate()
            .filter(|(index, val)| **val != Scalar::zero())
            .last()
            .map(|(index, _)| index)
            .unwrap_or(usize::MAX)
    }

    pub fn evaluate(&self, x: &Scalar) -> Scalar {
        let mut value = Scalar::zero();
        for term in self.terms.iter().rev() {
            value *= x;
            value += term;
        }
        value
    }
}

impl<const TERMS: usize> Mul for Polynomial<TERMS> {
    type Output = Polynomial<TERMS>;
    fn mul(self, rhs: Self) -> Self::Output {
        let deg_1 = self.degree();
        let deg_2 = rhs.degree();
        let new_degree = deg_1 + deg_2;
        assert!(new_degree < TERMS);
        let mut polynomial = Polynomial {
            terms: [Scalar::zero(); TERMS],
        };
        for i in 0..=new_degree {
            for lhs_index in i.saturating_sub(deg_2)..=deg_1.min(i) {
                let rhs_index = i - lhs_index;
                polynomial.terms[i] += self.terms[lhs_index] * rhs.terms[rhs_index];
            }
        }
        polynomial
    }
}

impl<const TERMS: usize> Mul<Scalar> for Polynomial<TERMS> {
    type Output = Self;

    fn mul(mut self, rhs: Scalar) -> Self::Output {
        for coef in self.terms.iter_mut() {
            *coef *= rhs;
        }
        self
    }
}

impl<const TERMS: usize> MulAssign<Polynomial<TERMS>> for Polynomial<TERMS> {
    fn mul_assign(&mut self, rhs: Polynomial<TERMS>) {
        *self = *self * rhs;
    }
}

impl<const TERMS: usize> MulAssign<Scalar> for Polynomial<TERMS> {
    fn mul_assign(&mut self, rhs: Scalar) {
        *self = *self * rhs;
    }
}

impl<const TERMS: usize> Mul<Polynomial<TERMS>> for Scalar {
    type Output = Polynomial<TERMS>;
    fn mul(self, rhs: Polynomial<TERMS>) -> Self::Output {
        rhs * self
    }
}

impl<const TERMS: usize> Add for Polynomial<TERMS> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        for (coef, other) in self.terms.iter_mut().zip(rhs.terms) {
            *coef += other;
        }
        self
    }
}

impl<const TERMS: usize> Sum for Polynomial<TERMS> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Polynomial::zero();
        for i in iter {
            res = res + i;
        }
        res
    }
}

pub mod lagrange {
    use bls12_381::Scalar;

    use super::Polynomial;

    /// Get an array of polynomials such that the i'th polynomial is 1 at
    /// evaluation point i and 0 on the other evaluation points.
    pub fn basis_polynomials<const N: usize>(evaluation_points: [Scalar; N]) -> [Polynomial<N>; N] {
        let mut result = [Polynomial::identity(); N];
        for i in 0..N {
            let numerator_polynomial = &mut result[i];
            let mut denominator = Scalar::one();
            for j in 0..N {
                if i != j {
                    *numerator_polynomial *= Polynomial::x_minus_const(evaluation_points[j]);
                    denominator *= evaluation_points[i] - evaluation_points[j];
                }
            }
            *numerator_polynomial *= denominator.invert().unwrap();
        }
        result
    }

    pub fn interpolate<const N: usize>(
        basis_polynomials: [Polynomial<N>; N],
        evaluations: [Scalar; N],
    ) -> Polynomial<N> {
        basis_polynomials
            .iter()
            .zip(evaluations.iter())
            .map(|(poly, val)| *val * *poly)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use bls12_381::Scalar;

    use crate::{
        polynomial::{lagrange, Polynomial},
        roots_of_unity,
    };

    #[test]
    fn test_evaluate() {
        // f(x) = x - 1
        let poly = Polynomial::<2>::x_minus_const(Scalar::one());
        assert_eq!(poly.evaluate(&Scalar::zero()), -Scalar::one());
        assert_eq!(poly.evaluate(&Scalar::one()), Scalar::zero());
    }

    #[test]
    fn test_lagrange_basis_polynomials() {
        const N: usize = 2 << 4;
        let roots_of_unity = roots_of_unity::roots_of_unity::<N>();
        let lagrange_basis_polynomials = lagrange::basis_polynomials(roots_of_unity);

        // Test whether they are correct
        for (i, poly) in lagrange_basis_polynomials.iter().enumerate() {
            for (j, evaluation_point) in roots_of_unity.iter().enumerate() {
                let val = poly.evaluate(evaluation_point);
                if i == j {
                    assert_eq!(val, Scalar::one());
                } else {
                    assert_eq!(val, Scalar::zero());
                }
            }
        }
    }
}
