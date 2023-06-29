use bls12_381::Scalar;

use crate::{
    polynomial::{lagrange, Polynomial},
    roots_of_unity::{self, roots_of_unity},
};

struct ConstraintMatrix<
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
    const CONSTRAINTS: usize,
> {
    public_entries: [[Scalar; CONSTRAINTS]; PUBLIC_WITNESS],
    private_entries: [[Scalar; CONSTRAINTS]; PRIVATE_WITNESS],
}

pub struct R1CS<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
{
    a: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    b: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    c: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
}

pub struct QAP<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const POLY_TERMS: usize> {
    pub public_v: [Polynomial<POLY_TERMS>; PUBLIC_WITNESS],
    pub private_v: [Polynomial<POLY_TERMS>; PRIVATE_WITNESS],
    pub public_w: [Polynomial<POLY_TERMS>; PUBLIC_WITNESS],
    pub private_w: [Polynomial<POLY_TERMS>; PRIVATE_WITNESS],
    pub public_y: [Polynomial<POLY_TERMS>; PUBLIC_WITNESS],
    pub private_y: [Polynomial<POLY_TERMS>; PRIVATE_WITNESS],
}

impl<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
    From<R1CS<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>>
    for QAP<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>
{
    fn from(r1cs: R1CS<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>) -> Self {
        let mut qap = QAP {
            public_v: [Polynomial::zero(); PUBLIC_WITNESS],
            public_w: [Polynomial::zero(); PUBLIC_WITNESS],
            public_y: [Polynomial::zero(); PUBLIC_WITNESS],
            private_v: [Polynomial::zero(); PRIVATE_WITNESS],
            private_w: [Polynomial::zero(); PRIVATE_WITNESS],
            private_y: [Polynomial::zero(); PRIVATE_WITNESS],
        };

        let roots_of_unity = roots_of_unity::<CONSTRAINTS>();
        let basis_polynomials = lagrange::basis_polynomials(roots_of_unity);

        for i in 0..PUBLIC_WITNESS {
            qap.public_v[i] = lagrange::interpolate(basis_polynomials, r1cs.a.public_entries[i]);
            qap.public_w[i] = lagrange::interpolate(basis_polynomials, r1cs.b.public_entries[i]);
            qap.public_y[i] = lagrange::interpolate(basis_polynomials, r1cs.c.public_entries[i]);
        }

        for i in 0..PRIVATE_WITNESS {
            qap.private_v[i] = lagrange::interpolate(basis_polynomials, r1cs.a.private_entries[i]);
            qap.private_w[i] = lagrange::interpolate(basis_polynomials, r1cs.b.private_entries[i]);
            qap.private_y[i] = lagrange::interpolate(basis_polynomials, r1cs.c.private_entries[i]);
        }

        qap
    }
}
