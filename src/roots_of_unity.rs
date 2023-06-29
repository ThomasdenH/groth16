use bls12_381::Scalar;
use ff::PrimeField;

/// Get the n'th root of unity. n should be a power of two.
pub fn primitive_root_of_unity(mut n: usize) -> Scalar {
    assert!(n.is_power_of_two());
    assert!(n <= (1 << 32));
    let mut root_of_unity = Scalar::ROOT_OF_UNITY;
    while n < (1 << 32) {
        n <<= 1;
        root_of_unity = root_of_unity.square();
    }
    root_of_unity
}

/// Get all 2^n powers of w, the 2^n'th root of unity.
pub fn roots_of_unity<const N: usize>() -> [Scalar; N] {
    let root_of_unity = primitive_root_of_unity(N);
    let mut all_roots = [Scalar::one(); N];
    all_roots[0] = root_of_unity;
    for i in 1..N {
        all_roots[i] = all_roots[i - 1] * root_of_unity;
    }
    all_roots
}

#[cfg(test)]
mod tests {
    use bls12_381::Scalar;

    use super::primitive_root_of_unity;

    #[test]
    fn test_primitive_root_of_unity() {
        for pow in 1usize..=32 {
            let two_pow = 1 << pow;
            let mut root_of_unity = primitive_root_of_unity(two_pow);
            for _ in 0..pow {
                assert_ne!(root_of_unity, Scalar::one());
                root_of_unity = root_of_unity.square();
            }
            assert_eq!(root_of_unity, Scalar::one());
        }
    }
}
