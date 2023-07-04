use std::{
    fmt::Debug,
    iter::{Sum, self},
    ops::{Add, Div, Mul, MulAssign, Sub},
};

use bls12_381::Scalar;
use pairing::PairingCurveAffine;

pub const ZERO: Polynomial = Polynomial { terms: Vec::new() };

#[derive(Clone, Eq, PartialEq, Default)]
pub struct Polynomial {
    terms: Vec<Scalar>,
}

impl Debug for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, coef) in self.terms.iter().enumerate().rev() {
            match i {
                0 => {
                    write!(f, "{coef}")?;
                }
                1 => {
                    write!(f, "{coef} * x + ")?;
                }
                i => {
                    write!(f, "{coef} * x^{i} + ")?;
                }
            }
        }
        Ok(())
    }
}

impl Polynomial {
    /// The polynomial t(x) = (x - w_1)(x - w_2)...
    /// for all roots of unity. The polynomial with 2^k + 1 terms
    /// will have its highest degree at 2^k as expected.
    pub fn t_x(degree: usize) -> Self {
        Polynomial {
            terms: iter::once(-Scalar::one())
                .chain(iter::repeat(Scalar::zero()).take(degree - 1))
                .chain(iter::once(Scalar::one()))
                .collect()
        }
    }

    pub const fn zero() -> Self {
        Polynomial {
            terms: Vec::new()
        }
    }

    /// The identity polynomial, equal to constant 1.
    /// 
    /// # Example
    /// ```rust
    /// assert_eq!(Polynomial::identity().evaluate(Scalar::one()), Scalar::one());
    /// assert_eq!(Polynomial::identity() * Polynomial::identity(), Polynomial::identity());
    /// ```
    pub fn identity() -> Self {
        Polynomial { terms: vec![Scalar::one()] }
    }

    /// Get the polynomial corresponding to (x - `const_term`).
    /// 
    /// # Example
    /// ```rust
    /// let poly = Polynomial::x_minus_const(5.into());
    /// assert_eq!(poly.evaluate(5.into()), Scalar::zero());
    /// ```
    pub fn x_minus_const(const_term: Scalar) -> Self {
        Polynomial { terms: vec![-const_term, Scalar::one()] }
    }

    /// Get the degree of this polynomial. Returns usize::MAX for zero-polynomials
    pub fn degree(&self) -> usize {
        self.terms.len().wrapping_sub(1)
    }

    pub fn evaluate(&self, x: &Scalar) -> Scalar {
        self.terms.iter()
            .rev()
            .fold(Scalar::zero(), |value, term| value * x + term)
    }

    pub fn evaluate_secret_powers<G, H>(&self, powers: &[G]) -> G
    where
        G: PairingCurveAffine + Mul<Scalar, Output = H>,
        H: Into<G> + Sum,
    {
        assert!(powers.len() >= self.terms.len());
        self.terms
            .iter()
            .zip(powers.iter())
            .map(|(coef, pow)| (*pow * *coef))
            .sum::<H>()
            .into()
    }

    /// Trim highest powers of x with zero coefficients
    fn trim(&mut self) {
        let to_truncate = self.terms.iter().rev().take_while(|coef| **coef == Scalar::zero()).count();
        self.terms.truncate(self.terms.len() - to_truncate);
    }

    fn highest_power_term(&self) -> Option<Scalar> {
        self.terms.last().copied()
    }

    fn is_trimmed(&self) -> bool {
        self.terms.last() != Some(&Scalar::zero())
    }

    fn is_zero(&self) -> bool {
        debug_assert!(self.is_trimmed());
        self.terms.is_empty()
    }
}

impl Mul for Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: Self) -> Self::Output {
        assert!(self.is_trimmed());
        assert!(rhs.is_trimmed());

        if self.is_zero() || rhs.is_zero() {
            return Polynomial::zero();
        }

        let degree_1 = self.degree();
        let degree_2 = rhs.degree();
        let new_degree = degree_1 + degree_2;
        let mut polynomial = Polynomial {
            terms: vec![Scalar::zero(); new_degree + 1]
        };
        for i in 0..=new_degree {
            for lhs_index in i.saturating_sub(degree_2)..=degree_1.min(i) {
                let rhs_index = i - lhs_index;
                polynomial.terms[i] += self.terms[lhs_index] * rhs.terms[rhs_index];
            }
        }
        polynomial.trim();
        polynomial
    }
}

impl Mul<Scalar> for Polynomial {
    type Output = Self;

    fn mul(mut self, rhs: Scalar) -> Self::Output {
        for coef in self.terms.iter_mut() {
            *coef *= rhs;
        }
        self.trim();
        self
    }
}

impl MulAssign<Polynomial> for Polynomial {
    fn mul_assign(&mut self, rhs: Polynomial) {
        *self = self.clone() * rhs;
    }
}

impl MulAssign<Scalar> for Polynomial {
    fn mul_assign(&mut self, rhs: Scalar) {
        *self = self.clone() * rhs;
    }
}

impl Mul<Polynomial> for Scalar {
    type Output = Polynomial;
    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs * self
    }
}

impl Add for Polynomial {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        for (coef, other) in self.terms.iter_mut().zip(rhs.terms) {
            *coef += other;
        }
        self.trim();
        self
    }
}

impl Sub for Polynomial {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        for (coef, other) in self.terms.iter_mut().zip(rhs.terms) {
            *coef -= other;
        }
        self.trim();
        self
    }
}

impl Sum for Polynomial {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Polynomial::zero();
        for i in iter {
            res = res + i;
        }
        res.trim();
        res
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct DivResult {
    pub quotient: Polynomial,
    pub remainder: Polynomial,
}

impl Div for &Polynomial {
    type Output = DivResult;
    fn div(self, rhs: Self) -> Self::Output {
        assert!(rhs != &Polynomial::zero());
        assert!(self.is_trimmed());
        assert!(rhs.is_trimmed());

        // At the end of this function,
        // self = rhs * quotient + remainder
        let rhs_degree = rhs.degree();
        let resulting_degree = self.degree() - rhs_degree;
        let mut remainder = self.clone();
        let mut quotient = Polynomial { terms: vec![Scalar::zero(); resulting_degree + 1]};

        while remainder != Polynomial::zero() && remainder.degree() >= rhs.degree() {
            // Divide leading coefficients
            let quotient_x_power = remainder.degree() - rhs_degree;
            quotient.terms[quotient_x_power] = remainder.highest_power_term().unwrap() * rhs.highest_power_term().unwrap().invert().unwrap();
            remainder = self.clone() - quotient.clone() * rhs.clone();
        }
        DivResult {
            quotient,
            remainder,
        }
    }
}

pub mod lagrange {
    use bls12_381::Scalar;

    use super::{Polynomial, ZERO};

    /// Get an array of polynomials such that the i'th polynomial is 1 at
    /// evaluation point i and 0 on the other evaluation points.
    pub fn basis_polynomials<const N: usize>(evaluation_points: [Scalar; N]) -> [Polynomial; N] {
        let mut result = [ZERO; N];
        for poly in result.iter_mut() {
            *poly = Polynomial::identity();
        }
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
        basis_polynomials: &[Polynomial; N],
        evaluations: &[Scalar; N],
    ) -> Polynomial {
        basis_polynomials
            .iter()
            .zip(evaluations.iter())
            .map(|(poly, val)| *val * poly.clone())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use bls12_381::Scalar;

    use crate::{
        polynomial::{lagrange, DivResult, Polynomial},
        roots_of_unity,
    };

    #[test]
    fn test_evaluate() {
        // f(x) = x - 1
        let poly = Polynomial::x_minus_const(Scalar::one());
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

    #[test]
    fn test_mul_by_zero() {
        let poly_1 = Polynomial::zero();
        let poly_2 = Polynomial { terms: vec![0.into(), 3.into(), 5.into()]};
        assert_eq!(poly_1.clone() * poly_2.clone(), Polynomial::zero());
        assert_eq!(poly_2.clone() * poly_1.clone(), Polynomial::zero());
    }

    #[test]
    fn test_div() {
        let poly_1 = Polynomial {
            terms: vec![1.into(), 0.into(), 1.into()],
        };
        let poly_2 = Polynomial {
            terms: vec![1.into()],
        };
        assert_eq!(
            &poly_1 / &poly_2,
            DivResult {
                quotient: poly_1,
                remainder: Polynomial::zero()
            }
        );
    }

    #[test]
    fn test_div_2() {
        let poly_1 = Polynomial {
            terms: vec![-(Scalar::from(3)), 10.into(), -(Scalar::from(5)), 3.into()],
        };
        let poly_2 = Polynomial {
            terms: vec![1.into(), 3.into()],
        };
        assert_eq!(
            &poly_1 / &poly_2,
            DivResult {
                quotient: Polynomial {
                    terms: vec![4.into(), -Scalar::from(2), 1.into()]
                },
                remainder: Polynomial {
                    terms: vec![-Scalar::from(7)]
                }
            }
        )
    }
}
