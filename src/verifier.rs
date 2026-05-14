use ark_bn254::Fr;
use ark_ff::{BigInteger, FftField, Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::{
    core::{hash_leaf, hash_pair, integer_from_hash, random_fr_from_hash},
    data::{ProverValues, QueryProof},
};

fn verify_merke_path(leaf: Fr, path: &Vec<[u8; 32]>, commited_root: [u8; 32], index: usize) {
    let current_hash = compute_merkle_path_root(leaf, path, index);

    assert_eq!(
        current_hash, commited_root,
        "Merkle path verification failed: expected root does not match computed root"
    );
}

fn compute_merkle_path_root(leaf: Fr, path: &Vec<[u8; 32]>, index: usize) -> [u8; 32] {
    let mut current_hash = hash_leaf(&leaf);
    let mut current_index = index;
    for sibling_hash in path {
        if current_index % 2 == 0 {
            current_hash = hash_pair(&current_hash, sibling_hash);
        } else {
            current_hash = hash_pair(sibling_hash, &current_hash);
        }
        current_index /= 2;
    }
    current_hash
}

fn verify_trace_to_cp_consistency(
    query: Fr,
    f_x: Fr,
    f_gx: Fr,
    f_g2x: Fr,
    cp_x: Fr,
    trace_length: usize,
    _blowup_factor: usize,
    target_value: Fr,
) {
    let domain = GeneralEvaluationDomain::<Fr>::new(trace_length)
        .unwrap()
        .elements()
        .collect::<Vec<Fr>>();
    let g = domain[1];

    let p0_at_x = (f_x - Fr::from(1u64)) / (query - Fr::from(1u64));
    let p1_at_x = (f_x - target_value) / (query - g.pow([trace_length as u64 - 1u64]));
    let p2_at_x = (f_g2x - f_gx.square() - f_x.square())
        / (query.pow([trace_length as u64]) - Fr::from(1u64))
        * (query - g.pow([trace_length as u64 - 1u64]))
        * (query - g.pow([trace_length as u64 - 2u64]));

    let mut alpha_seed: Vec<u8> = Vec::new();
    // use the target and generator of the domain as part of the seed
    alpha_seed.extend_from_slice(target_value.into_bigint().to_bytes_le().as_slice());
    alpha_seed.extend_from_slice(g.into_bigint().to_bytes_le().as_slice());
    let alpha_0 = random_fr_from_hash(&alpha_seed);
    alpha_seed.extend_from_slice(alpha_0.into_bigint().to_bytes_le().as_slice());
    let alpha_1 = random_fr_from_hash(&alpha_seed);
    alpha_seed.extend_from_slice(alpha_1.into_bigint().to_bytes_le().as_slice());
    let alpha_2 = random_fr_from_hash(&alpha_seed);

    let cp_at_x = &(&(&p0_at_x * &alpha_0) + &(&p1_at_x * &alpha_1)) + &(&p2_at_x * &alpha_2);

    assert_eq!(
        cp_at_x, cp_x,
        "CP consistency check failed: expected CP evaluation does not match computed CP evaluation"
    );
}

fn verify_fri_layer(
    previous_cp_at_x: &Fr,
    previous_cp_at_minux_x: &Fr,
    exptected_next_cp: &Fr,
    beta: &Fr,
    query: &Fr,
) {
    let two_inv = Fr::from(2u64).inverse().unwrap();
    let g_at_query = (previous_cp_at_x + previous_cp_at_minux_x) * two_inv;
    let h_at_query =
        (previous_cp_at_x - previous_cp_at_minux_x) * two_inv * query.inverse().unwrap();
    let next_cp_at_query = g_at_query + beta * &h_at_query;

    assert_eq!(
        next_cp_at_query, *exptected_next_cp,
        "FRI layer verification failed: expected next CP evaluation does not match computed next CP evaluation"
    );
}

fn verify_final_fri_layer(final_cp_x: Fr, final_cp_minus_x: Fr) {
    assert_eq!(
        final_cp_x, final_cp_minus_x,
        "Final FRI layer verification failed: final layer must be constant"
    );
}

fn verify_query(
    q: usize,

    extended_trace_commitment: &[u8; 32],
    composite_polynomial_commitment: &[u8; 32],
    fri_commitments: &Vec<[u8; 32]>,
    query_proof: &QueryProof,
    size: usize,
    blowup_factor: usize,
    target_value: Fr,
) {
    let extended_domain = GeneralEvaluationDomain::<Fr>::new(size * blowup_factor)
        .unwrap()
        .elements()
        .collect::<Vec<Fr>>();
    let query = Fr::GENERATOR * extended_domain[q];
    verify_merke_path(
        query_proof.f_x.value,
        &query_proof.f_x.path,
        *extended_trace_commitment,
        q,
    );
    let gq = (q + blowup_factor) % (size * blowup_factor);
    verify_merke_path(
        query_proof.f_gx.value,
        &query_proof.f_gx.path,
        *extended_trace_commitment,
        gq,
    );
    let ggq = (gq + blowup_factor) % (size * blowup_factor);
    verify_merke_path(
        query_proof.f_g2x.value,
        &query_proof.f_g2x.path,
        *extended_trace_commitment,
        ggq,
    );

    verify_trace_to_cp_consistency(
        query,
        query_proof.f_x.value,
        query_proof.f_gx.value,
        query_proof.f_g2x.value,
        query_proof.cp_pairs[0].cp_x.value,
        size,
        blowup_factor,
        target_value,
    );

    let mut beta_seed: Vec<u8> = Vec::new();
    beta_seed.extend_from_slice(extended_trace_commitment);
    beta_seed.extend_from_slice(composite_polynomial_commitment);
    let domain_size: usize = size * blowup_factor;
    let mut layer_index = q;
    let mut transition_query = query;
    for (i, cp_pair) in query_proof.cp_pairs.iter().enumerate() {
        let layer_size = domain_size >> i;
        if i != query_proof.cp_pairs.len() - 1 {
            let layer_commitment = if i == 0 {
                *composite_polynomial_commitment
            } else {
                fri_commitments[i - 1]
            };
            let negative_index = (layer_size / 2 + layer_index) % layer_size;
            verify_merke_path(
                cp_pair.cp_x.value,
                &cp_pair.cp_x.path,
                layer_commitment,
                layer_index,
            );
            verify_merke_path(
                cp_pair.cp_minus_x.value,
                &cp_pair.cp_minus_x.path,
                layer_commitment,
                negative_index,
            );
        }

        let beta = random_fr_from_hash(&beta_seed);
        if i != 0 {
            verify_fri_layer(
                &query_proof.cp_pairs[i - 1].cp_x.value,
                &query_proof.cp_pairs[i - 1].cp_minus_x.value,
                &cp_pair.cp_x.value,
                &beta,
                &transition_query,
            );
            transition_query = transition_query.square();
        }

        if i > 0 && i - 1 < fri_commitments.len() {
            beta_seed.extend_from_slice(&fri_commitments[i - 1]);
        }

        if i + 1 < query_proof.cp_pairs.len() {
            let next_layer_size = domain_size >> (i + 1);
            layer_index %= next_layer_size;
        }
    }
    let final_cp_pair = query_proof.cp_pairs.last().unwrap();
    verify_final_fri_layer(final_cp_pair.cp_x.value, final_cp_pair.cp_minus_x.value);
}

pub fn verify(
    prover_values: ProverValues,
    size: usize,
    blowup_factor: usize,
    _num_queries: usize,
    target_value: Fr,
) {
    let mut fiat_shamir_seed: Vec<u8> = Vec::new();
    fiat_shamir_seed.extend_from_slice(&prover_values.extended_trace_commitment);
    fiat_shamir_seed.extend_from_slice(&prover_values.composite_polynomial_commitment);
    for commitment in &prover_values.fri_commitments {
        fiat_shamir_seed.extend_from_slice(commitment);
    }
    for query_proof in prover_values.queries_proofs.iter() {
        let q = integer_from_hash(&fiat_shamir_seed, size * blowup_factor);
        verify_query(
            q,
            &prover_values.extended_trace_commitment,
            &prover_values.composite_polynomial_commitment,
            &prover_values.fri_commitments,
            query_proof,
            size,
            blowup_factor,
            target_value,
        );

        for node in query_proof.f_x.path.iter() {
            fiat_shamir_seed.extend_from_slice(node);
        }
        for node in query_proof.f_gx.path.iter() {
            fiat_shamir_seed.extend_from_slice(node);
        }
        for node in query_proof.f_g2x.path.iter() {
            fiat_shamir_seed.extend_from_slice(node);
        }
        for cp_pair in query_proof.cp_pairs.iter() {
            for node in cp_pair.cp_x.path.iter() {
                fiat_shamir_seed.extend_from_slice(node);
            }
            for node in cp_pair.cp_minus_x.path.iter() {
                fiat_shamir_seed.extend_from_slice(node);
            }
        }
    }
}
