use ark_bn254::Fr;

pub struct ValueAndPath {
    pub value: Fr,
    pub path: Vec<[u8; 32]>,
}

pub struct CpPair {
    pub cp_x: ValueAndPath,
    pub cp_minus_x: ValueAndPath,
}

pub struct QueryProof {
    pub f_x: ValueAndPath,
    pub f_gx: ValueAndPath,
    pub f_g2x: ValueAndPath,
    pub cp_pairs: Vec<CpPair>,
}

pub struct ProverValues {
    pub extended_trace_commitment: [u8; 32],
    pub composite_polynomial_commitment: [u8; 32],
    pub fri_commitments: Vec<[u8; 32]>,
    pub queries_proofs: Vec<QueryProof>,
}
