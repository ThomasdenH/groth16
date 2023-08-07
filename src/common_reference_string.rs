use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use pairing::PairingCurveAffine;

use crate::{polynomial::Polynomial, r1cs::Qap};

/// Get an array with powers of x on a curve, i.e. obtain G, xG, x^2G, etc.
/// Uses the default generator.
fn powers_of_x<
    const N: usize,
    Curve: PairingCurveAffine
        + std::ops::Mul<bls12_381::Scalar>
        + std::convert::From<<Curve as std::ops::Mul<bls12_381::Scalar>>::Output>,
>(
    x: Scalar,
) -> [Curve; N] {
    powers_of_x_with_generator(x, Curve::generator())
}

/// Get an array with powers of x on a curve. Given x and G, this will return
/// `[G, xG, ...]`.
fn powers_of_x_with_generator<
    const N: usize,
    Curve: PairingCurveAffine
        + std::ops::Mul<bls12_381::Scalar>
        + std::convert::From<<Curve as std::ops::Mul<bls12_381::Scalar>>::Output>,
>(
    x: Scalar,
    generator: Curve,
) -> [Curve; N] {
    let mut powers = [generator; N];
    for i in 1..N {
        powers[i] = (powers[i - 1] * x).into();
    }
    powers
}

/// A common reference string for a particular computation.
pub struct CommonReferenceString<
    const N: usize, // Highest power of x
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
> {
    /// The CRS elements in the first curve group.
    g1: CommonReferenceStringG1<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>,
    /// The CRS elements in the second curve group.
    g2: CommonReferenceStringG2<N>,
}

struct CommonReferenceStringG1<
    const N: usize,
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
> {
    alpha: G1Affine,
    beta: G1Affine,
    delta: G1Affine,
    x_power: [G1Affine; N],
    public_contribs: [G1Affine; PUBLIC_WITNESS],
    private_contribs: [G1Affine; PRIVATE_WITNESS],
    x_power_t_x: [G1Affine; N_MINUS_ONE],
}

struct CommonReferenceStringG2<const N: usize> {
    beta: G2Affine,
    gamma: G2Affine,
    delta: G2Affine,
    x_power: [G2Affine; N],
}

impl<
        const N: usize, // Highest power of x
        const N_MINUS_ONE: usize,
        const PUBLIC_WITNESS: usize,
        const PRIVATE_WITNESS: usize,
    > CommonReferenceString<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>
{
    /// Build a new common reference string from the given trapdoor and
    /// Quadratic Arithmetic Program.
    pub fn new<const N_PLUS_ONE: usize>(
        alpha: Scalar,
        beta: Scalar,
        gamma: Scalar,
        delta: Scalar,
        x: Scalar,
        qap: &Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>,
    ) -> Self {
        assert_eq!(N_PLUS_ONE, N + 1);
        assert_eq!(N_MINUS_ONE, N - 1);
        // t(x)
        let t_x = Polynomial::t_x(N).evaluate(&x);
        // t(x) / delta
        let t_x_over_delta = t_x * delta.invert().unwrap() * G1Affine::generator();
        // x^i t(x) / delta
        let x_power_t_x = powers_of_x_with_generator(x, t_x_over_delta.into());

        let mut public_contribs = [G1Affine::generator(); PUBLIC_WITNESS];
        for (i, public_contrib) in public_contribs.iter_mut().enumerate() {
            *public_contrib = ((beta * qap.public_v[i].evaluate(&x)
                + alpha * qap.public_w[i].evaluate(&x)
                + qap.public_y[i].evaluate(&x))
                * gamma.invert().unwrap()
                * G1Affine::generator())
            .into();
        }
        let mut private_contribs = [G1Affine::generator(); PRIVATE_WITNESS];
        for (i, private_contrib) in private_contribs.iter_mut().enumerate() {
            *private_contrib = ((beta * qap.private_v[i].evaluate(&x)
                + alpha * qap.private_w[i].evaluate(&x)
                + qap.private_y[i].evaluate(&x))
                * delta.invert().unwrap()
                * G1Affine::generator())
            .into();
        }

        CommonReferenceString {
            g1: CommonReferenceStringG1 {
                alpha: (G1Affine::generator() * alpha).into(),
                beta: (G1Affine::generator() * beta).into(),
                delta: (G1Affine::generator() * delta).into(),
                x_power: powers_of_x(x),
                x_power_t_x,
                public_contribs,
                private_contribs,
            },
            g2: CommonReferenceStringG2 {
                beta: (G2Affine::generator() * beta).into(),
                delta: (G2Affine::generator() * delta).into(),
                gamma: (G2Affine::generator() * gamma).into(),
                x_power: powers_of_x(x),
            },
        }
    }
}

pub struct Proof {
    a: G1Affine,
    b: G2Affine,
    c: G1Affine,
}

impl Proof {
    pub fn verify<
        const N: usize, // Highest power of x
        const N_MINUS_ONE: usize,
        const PUBLIC_WITNESS: usize,
        const PRIVATE_WITNESS: usize,
    >(
        &self,
        crs: &CommonReferenceString<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>,
        public_assignment: &[Scalar; PUBLIC_WITNESS],
    ) -> bool {
        self.a.pairing_with(&self.b)
        // Can be precomputed:
        == crs.g1.alpha.pairing_with(&crs.g2.beta)
        + G1Affine::from(public_assignment.iter()
            .zip(crs.g1.public_contribs.iter())
            .map(|(a, contrib)| a * contrib)
            .sum::<G1Projective>()).pairing_with(&crs.g2.gamma)
        + self.c.pairing_with(&crs.g2.delta)
    }
}

/// Build a proof using a common reference string, a QAP and an assignment.
pub fn prove<
    const N: usize, // Highest power of x
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
>(
    common_reference_string: &CommonReferenceString<
        N,
        N_MINUS_ONE,
        PUBLIC_WITNESS,
        PRIVATE_WITNESS,
    >,
    qap: Qap<PUBLIC_WITNESS, PRIVATE_WITNESS>,
    public_assignment: [Scalar; PUBLIC_WITNESS],
    private_assignment: [Scalar; PRIVATE_WITNESS],
) -> Proof {
    let r = Scalar::one();
    let s = Scalar::one();

    let v_x: Polynomial = public_assignment
        .iter()
        .zip(qap.public_v.iter())
        .map(|(a, v)| *a * v.clone())
        .sum::<Polynomial>()
        + private_assignment
            .iter()
            .zip(qap.private_v.iter())
            .map(|(a, v)| *a * v.clone())
            .sum();

    let w_x: Polynomial = public_assignment
        .iter()
        .zip(qap.public_w.iter())
        .map(|(a, w)| *a * w.clone())
        .sum::<Polynomial>()
        + private_assignment
            .iter()
            .zip(qap.private_w.iter())
            .map(|(a, w)| *a * w.clone())
            .sum();

    let y_x: Polynomial = public_assignment
        .iter()
        .zip(qap.public_y.iter())
        .map(|(a, y)| *a * y.clone())
        .sum::<Polynomial>()
        + private_assignment
            .iter()
            .zip(qap.private_y.iter())
            .map(|(a, y)| *a * y.clone())
            .sum();

    let p_x = v_x.clone() * w_x.clone() - y_x;
    let q_x = &p_x / &Polynomial::t_x(N);
    assert!(q_x.remainder == Polynomial::zero());

    let a = (common_reference_string.g1.alpha
        + G1Projective::from(v_x.evaluate_secret_powers(&common_reference_string.g1.x_power))
        + r * common_reference_string.g1.delta)
        .into();
    let b = (common_reference_string.g2.beta
        + G2Projective::from(w_x.evaluate_secret_powers(&common_reference_string.g2.x_power))
        + s * common_reference_string.g2.delta)
        .into();
    let b_g1: G1Affine = (common_reference_string.g1.beta
        + G1Projective::from(w_x.evaluate_secret_powers(&common_reference_string.g1.x_power))
        + s * common_reference_string.g1.delta)
        .into();
    let c = G1Affine::from(
        private_assignment
            .iter()
            .zip(common_reference_string.g1.private_contribs.iter())
            .map(|(a, contrib)| *a * *contrib)
            .sum::<G1Projective>()
            + q_x
                .quotient
                .evaluate_secret_powers(&common_reference_string.g1.x_power_t_x)
            + a * s
            + r * b_g1
            - r * s * common_reference_string.g1.delta,
    );

    Proof { a, b, c }
}

#[cfg(test)]
mod tests {
    use crate::r1cs::{self, Qap};

    use super::{prove, CommonReferenceString};

    #[test]
    fn test_prove() {
        let qap: Qap<1, 1> = r1cs::tests::R1CS.into();
        let crs = CommonReferenceString::<2, 1, 1, 1>::new::<3>(
            123.into(),
            2346.into(),
            2136.into(),
            15.into(),
            4756.into(),
            &qap,
        );
        let proof = prove(
            &crs,
            qap,
            r1cs::tests::R1CS_ASSIGNMENT_PUBLIC,
            r1cs::tests::R1CS_ASSIGNMENT_PRIVATE,
        );
        assert!(proof.verify(&crs, &r1cs::tests::R1CS_ASSIGNMENT_PUBLIC));
    }
}
