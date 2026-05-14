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

pub fn hash_leaf(value: &Fr) -> [u8; 32] {
    let bytes = fr_to_bytes(value);
    Sha256::digest(bytes)
        .try_into()
        .expect("Hash output should be 32 bytes")
}

pub fn compute_merkle_root(leaves_values: &[Fr]) -> (Vec<Vec<[u8; 32]>>, [u8; 32]) {
    assert!(
        !leaves_values.is_empty(),
        "merkle tree requires at least one leaf"
    );

    let mut current_level: Vec<[u8; 32]> = leaves_values.iter().map(hash_leaf).collect();
    let mut merkle_tree: Vec<Vec<[u8; 32]>> = Vec::new();

    while current_level.len() > 1 {
        merkle_tree.push(current_level.clone());
        current_level = current_level
            .chunks(2)
            .map(|chunk| {
                let left = chunk[0];
                let right = if chunk.len() > 1 { chunk[1] } else { [0u8; 32] };
                hash_pair(&left, &right)
            })
            .collect();
    }
    merkle_tree.push(current_level.clone());
    (merkle_tree, current_level[0])
}

pub fn integer_from_hash(hash: &[u8], modulus: usize) -> usize {
    if modulus == 0 {
        return 0;
    }
    let mut result: usize = 0;
    for byte in hash {
        result = (result.wrapping_mul(256).wrapping_add(*byte as usize)) % modulus;
    }
    result
}

pub fn random_fr_from_hash(input: &[u8]) -> Fr {
    let hash = sha2::Sha256::digest(input);
    let mut bytes = [0u8; 32];
    bytes.copy_from_slice(&hash);
    Fr::from_le_bytes_mod_order(&bytes)
}
