use ark_bn254::Fr;

pub struct ProverValues {
    pub extended_trace_commitment: [u8; 32],
    pub composite_polynomial_commitment: [u8; 32],
    pub fri_commitments: Vec<[u8; 32]>,
}
