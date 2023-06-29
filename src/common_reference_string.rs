use bls12_381::{G1Affine, G2Affine, Scalar};
use pairing::PairingCurveAffine;

use crate::{polynomial::Polynomial, r1cs::QAP};

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

struct CommonReferenceString<
    const N: usize, // Highest power of x
    const N_MINUS_ONE: usize,
    const PUBLIC_WITNESS: usize,
    const PRIVATE_WITNESS: usize,
> {
    g1: CommonReferenceStringG1<N, N_MINUS_ONE, PUBLIC_WITNESS, PRIVATE_WITNESS>,
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
    fn new<const N_PLUS_ONE: usize>(
        alpha: Scalar,
        beta: Scalar,
        gamma: Scalar,
        delta: Scalar,
        x: Scalar,
        qap: QAP<PUBLIC_WITNESS, PRIVATE_WITNESS, N>,
    ) -> Self {
        assert_eq!(N_PLUS_ONE, N + 1);
        assert_eq!(N_MINUS_ONE, N - 1);
        // t(x)
        let t_x = Polynomial::<N_PLUS_ONE>::t_x(N).evaluate(&x);
        // t(x) / delta
        let t_x_over_delta = t_x * delta.invert().unwrap() * G1Affine::generator();
        // x^i t(x) / delta
        let x_power_t_x = powers_of_x_with_generator(x, t_x_over_delta.into());

        let mut public_contribs = [G1Affine::generator(); PUBLIC_WITNESS];
        for i in 0..PUBLIC_WITNESS {
            public_contribs[i] = ((beta * qap.public_v[i].evaluate(&x)
                + alpha * qap.public_w[i].evaluate(&x)
                + qap.public_y[i].evaluate(&x))
                * gamma.invert().unwrap()
                * G1Affine::generator())
            .into();
        }
        let mut private_contribs = [G1Affine::generator(); PRIVATE_WITNESS];
        for i in 0..PRIVATE_WITNESS {
            public_contribs[i] = ((beta * qap.private_v[i].evaluate(&x)
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
