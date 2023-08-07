use bls12_381::Scalar;

use crate::{
    polynomial::{self, lagrange, Polynomial},
    roots_of_unity::roots_of_unity,
};

pub struct ConstraintMatrix<
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
    const CONSTRAINTS: usize,
> {
    public_entries: [[Scalar; CONSTRAINTS]; PUBLIC_WITNESS],
    private_entries: [[Scalar; CONSTRAINTS]; PRIVATE_WITNESS],
}

pub struct R1cs<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
{
    a: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    b: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    c: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
}

impl<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
    R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>
{
    pub const fn new(
        a: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
        b: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
        c: ConstraintMatrix<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>,
    ) -> Self {
        R1cs { a, b, c }
    }

    pub fn is_satisfied_by(
        &self,
        public_witness: [Scalar; PUBLIC_WITNESS],
        private_witness: [Scalar; PRIVATE_WITNESS],
    ) -> bool {
        for constraint_index in 0..CONSTRAINTS {
            let mut a_sum = Scalar::zero();
            let mut b_sum = Scalar::zero();
            let mut c_sum = Scalar::zero();
            for (private_idx, private_witness_entry) in private_witness.iter().enumerate() {
                a_sum +=
                    private_witness_entry * self.a.private_entries[private_idx][constraint_index];
                b_sum +=
                    private_witness_entry * self.b.private_entries[private_idx][constraint_index];
                c_sum +=
                    private_witness_entry * self.c.private_entries[private_idx][constraint_index];
            }
            for (witness_idx, public_witness_entry) in public_witness.iter().enumerate() {
                a_sum +=
                    public_witness_entry * self.a.public_entries[witness_idx][constraint_index];
                b_sum +=
                    public_witness_entry * self.b.public_entries[witness_idx][constraint_index];
                c_sum +=
                    public_witness_entry * self.c.public_entries[witness_idx][constraint_index];
            }
            if a_sum * b_sum != c_sum {
                return false;
            }
        }
        true
    }
}

#[derive(Clone)]
pub struct Qap<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize> {
    pub public_v: [Polynomial; PUBLIC_WITNESS],
    pub private_v: [Polynomial; PRIVATE_WITNESS],
    pub public_w: [Polynomial; PUBLIC_WITNESS],
    pub private_w: [Polynomial; PRIVATE_WITNESS],
    pub public_y: [Polynomial; PUBLIC_WITNESS],
    pub private_y: [Polynomial; PRIVATE_WITNESS],
}

impl<const PUBLIC_WITNESS: usize, const PRIVATE_WITNESS: usize, const CONSTRAINTS: usize>
    From<R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>>
    for Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>
{
    fn from(r1cs: R1cs<PUBLIC_WITNESS, PRIVATE_WITNESS, CONSTRAINTS>) -> Self {
        let mut qap = Qap {
            public_v: [polynomial::ZERO; PUBLIC_WITNESS],
            public_w: [polynomial::ZERO; PUBLIC_WITNESS],
            public_y: [polynomial::ZERO; PUBLIC_WITNESS],
            private_v: [polynomial::ZERO; PRIVATE_WITNESS],
            private_w: [polynomial::ZERO; PRIVATE_WITNESS],
            private_y: [polynomial::ZERO; PRIVATE_WITNESS],
        };

        let roots_of_unity = roots_of_unity::<CONSTRAINTS>();
        let basis_polynomials = lagrange::basis_polynomials(roots_of_unity);

        for i in 0..PUBLIC_WITNESS {
            qap.public_v[i] = lagrange::interpolate(&basis_polynomials, &r1cs.a.public_entries[i]);
            qap.public_w[i] = lagrange::interpolate(&basis_polynomials, &r1cs.b.public_entries[i]);
            qap.public_y[i] = lagrange::interpolate(&basis_polynomials, &r1cs.c.public_entries[i]);
        }

        for i in 0..PRIVATE_WITNESS {
            qap.private_v[i] =
                lagrange::interpolate(&basis_polynomials, &r1cs.a.private_entries[i]);
            qap.private_w[i] =
                lagrange::interpolate(&basis_polynomials, &r1cs.b.private_entries[i]);
            qap.private_y[i] =
                lagrange::interpolate(&basis_polynomials, &r1cs.c.private_entries[i]);
        }

        qap
    }
}

#[cfg(test)]
pub mod tests {
    use bls12_381::Scalar;

    use super::{ConstraintMatrix, R1cs};

    pub const R1CS: R1cs<1, 1, 2> = R1cs::new(
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::one()]],
            public_entries: [[Scalar::zero(), Scalar::zero()]],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::one(), Scalar::zero()]],
            public_entries: [[Scalar::zero(), Scalar::one()]],
        },
        ConstraintMatrix {
            private_entries: [[Scalar::zero(), Scalar::one()]],
            public_entries: [[Scalar::one(), Scalar::zero()]],
        },
    );

    pub const R1CS_ASSIGNMENT_PRIVATE: [Scalar; 1] = [Scalar::one()];
    pub const R1CS_ASSIGNMENT_PUBLIC: [Scalar; 1] = [Scalar::one()];

    #[test]
    fn test_r1cs() {
        // Use witness of size 2:
        // Private: a
        // Public: b
        // a^2 = b
        // ab = a
        assert!(R1CS.is_satisfied_by(R1CS_ASSIGNMENT_PRIVATE, R1CS_ASSIGNMENT_PUBLIC));
    }
}
