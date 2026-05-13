use ark_bn254::Fr;
use ark_ff::{BigInteger, PrimeField};
use sha2::{Digest, Sha256};

pub fn fr_to_bytes(fr: &Fr) -> [u8; 32] {
    let bigint = fr.into_bigint();
    let bytes = bigint.to_bytes_le();
    bytes
        .try_into()
        .expect("Fr element should fit into 32 bytes")
}
pub fn hash_pair(left: &[u8; 32], right: &[u8; 32]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(left);
    hasher.update(right);
    hasher
        .finalize()
        .try_into()
        .expect("Hash output should be 32 bytes")
}

pub fn compute_merkle_root(leaves_values: &Vec<Fr>) -> [u8; 32] {
    let mut current_level: Vec<[u8; 32]> = leaves_values.iter().map(fr_to_bytes).collect();
    while current_level.len() > 1 {
        current_level = current_level
            .chunks(2)
            .map(|chunk| {
                let left = chunk[0];
                let right = if chunk.len() > 1 { chunk[1] } else { [0u8; 32] };
                hash_pair(&left, &right)
            })
            .collect();
    }
    current_level[0]
}
