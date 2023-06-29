use std::{iter, ops::Mul};

use bls12_381::Scalar;
use pairing::PairingCurveAffine;

/// Powers of some secret point `X`.
pub struct SecretPointPowers<const TERMS: usize, Point: PairingCurveAffine> {
    terms: [Point; TERMS],
}

impl<const TERMS: usize, Point: PairingCurveAffine, MulOutput> SecretPointPowers<TERMS, Point>
where
    MulOutput: Into<Point>,
    Point: Mul<Scalar, Output = MulOutput>,
{
    pub fn new(x: Scalar) -> Self {
        SecretPointPowers {
            terms: iter::successors(Some(Point::generator()), |succ| Some((*succ * x).into()))
                .take(TERMS)
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
        }
    }
}
