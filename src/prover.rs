use ark_bn254::Fr;
use ark_ff::{BigInteger, FftField, Field, PrimeField, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DenseOrSparsePolynomial,
    univariate::DensePolynomial,
};

use crate::core::integer_from_hash;
use crate::core::random_fr_from_hash;
use crate::core::{compute_merkle_root, next_power_of_two};
use crate::data::CpPair;
use crate::data::ProverValues;
use crate::data::QueryProof;
use crate::data::ValueAndPath;

fn calculate_fib_square_trace(a0: Fr, a1: Fr, length: usize) -> Vec<Fr> {
    let mut trace: Vec<Fr> = vec![a0, a1];

    for i in 2..length {
        let next = trace[i - 1].square() + trace[i - 2].square();
        // println!("fib squared at index {}: {} ", i, next);
        trace.push(next);
    }

    trace
}

fn calculate_trace_polynomial(
    trace: &Vec<Fr>,
) -> (DensePolynomial<Fr>, GeneralEvaluationDomain<Fr>) {
    let domain = GeneralEvaluationDomain::<Fr>::new(trace.len()).unwrap();

    let evaluations = Evaluations::from_vec_and_domain(trace.to_vec(), domain);

    let polynomial = evaluations.interpolate();

    (polynomial, domain)
}

fn lde_on_coset(
    polynomial: &DensePolynomial<Fr>,
    size: usize,
    blowup_factor: usize,
    offset: Fr, // the value w in w·h^i
) -> (Vec<Fr>, GeneralEvaluationDomain<Fr>) {
    let extended_domain = GeneralEvaluationDomain::<Fr>::new(size * blowup_factor).unwrap();

    let shifted_coeffs = polynomial
        .coeffs
        .iter()
        .enumerate()
        .map(|(i, &c)| c * offset.pow(&[i as u64]))
        .collect::<Vec<_>>();

    (extended_domain.fft(&shifted_coeffs), extended_domain)
}

fn calculate_composite_polynomial(
    f: &DensePolynomial<Fr>,
    domain: &GeneralEvaluationDomain<Fr>,
    extended_domain: &GeneralEvaluationDomain<Fr>,
    offset: Fr,
    target: Fr,
    target_index: usize,
    trace_length: usize,
) -> (Vec<Fr>, Vec<Vec<[u8; 32]>>, [u8; 32]) {
    // let f: DensePolynomial<Fr> = DensePolynomial::from_coefficients_vec(domain.ifft(trace));

    let g = domain.group_gen();
    // calculating p0 = (f(X) - f(1)) / (X - 1)
    let p0_numerator: DenseOrSparsePolynomial<Fr> = DenseOrSparsePolynomial::from(
        f - &DensePolynomial::from_coefficients_vec(vec![Fr::from(1u64)]),
    );
    let p0_denominator: DenseOrSparsePolynomial<Fr> =
        DenseOrSparsePolynomial::from(DensePolynomial::from_coefficients_vec(vec![
            -Fr::from(1u64),
            Fr::from(1u64),
        ]));

    let (p0, p0_remainder) = p0_numerator.divide_with_q_and_r(&p0_denominator).unwrap();
    assert!(p0_remainder.is_zero(), "p0 division must be exact");

    // calculating p1 = (f(X) - a_{traget_index}) / (X - g^{target_index})
    let p1_numerator: DenseOrSparsePolynomial<Fr> =
        DenseOrSparsePolynomial::from(f - &DensePolynomial::from_coefficients_vec(vec![target]));
    let p1_denominator: DenseOrSparsePolynomial<Fr> =
        DenseOrSparsePolynomial::from(DensePolynomial::from_coefficients_vec(vec![
            -g.pow(&[target_index as u64]),
            Fr::from(1u64),
        ]));
    let (p1, p1_remainder) = p1_numerator.divide_with_q_and_r(&p1_denominator).unwrap();
    assert!(p1_remainder.is_zero(), "p1 division must be exact");

    // calculating p3 = (f((g^2)X) - f(gX)^2 - f(X)^2) * (X - g^{target_index}) * (X - g^{target_index-1}) / (X^{{target_index+1}} - 1)
    let f_g2x = DensePolynomial::from_coefficients_vec(
        f.coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * g.pow(&[2 * i as u64]))
            .collect(),
    );
    let f_gx = DensePolynomial::from_coefficients_vec(
        f.coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * g.pow(&[i as u64]))
            .collect(),
    );
    let f_gx_square = &f_gx * &f_gx;
    let f_square = f * f;
    let tmp = &f_g2x - &f_gx_square;
    let transition_numerator = &tmp - &f_square;

    let x_minus_g_target_index = DensePolynomial::from_coefficients_vec(vec![
        -g.pow(&[target_index as u64]),
        Fr::from(1u64),
    ]);
    let x_minus_g_target_index_minus_one = DensePolynomial::from_coefficients_vec(vec![
        -g.pow(&[(target_index - 1) as u64]),
        Fr::from(1u64),
    ]);
    let mut exclusion_points_polynomial =
        DensePolynomial::from_coefficients_vec(vec![Fr::from(1u64)]);
    for i in (target_index + 1)..trace_length {
        let factor =
            DensePolynomial::from_coefficients_vec(vec![-g.pow(&[i as u64]), Fr::from(1u64)]);
        exclusion_points_polynomial = &exclusion_points_polynomial * &factor;
    }
    let p3_numerator = DenseOrSparsePolynomial::from(
        &(&(&transition_numerator * &x_minus_g_target_index) * &x_minus_g_target_index_minus_one)
            * &exclusion_points_polynomial,
    );

    let p3_denominator = DenseOrSparsePolynomial::from(domain.vanishing_polynomial());
    let (p3, p3_remainder) = p3_numerator.divide_with_q_and_r(&p3_denominator).unwrap();
    assert!(p3_remainder.is_zero(), "p3 division must be exact");

    // now calculating CP using p0, p1, p3 and random linear combination
    let mut alpha_seed: Vec<u8> = Vec::new();
    // use the target and generator of the domain as part of the seed
    alpha_seed.extend_from_slice(target.into_bigint().to_bytes_le().as_slice());
    alpha_seed.extend_from_slice(g.into_bigint().to_bytes_le().as_slice());
    let alpha_0 = random_fr_from_hash(&alpha_seed);
    alpha_seed.extend_from_slice(alpha_0.into_bigint().to_bytes_le().as_slice());
    let alpha_1 = random_fr_from_hash(&alpha_seed);
    alpha_seed.extend_from_slice(alpha_1.into_bigint().to_bytes_le().as_slice());
    let alpha_2 = random_fr_from_hash(&alpha_seed);

    let composite_polynomial = &(&(&p0 * alpha_0) + &(&p1 * alpha_1)) + &(&p3 * alpha_2);
    let shifted_coeffs = composite_polynomial
        .coeffs
        .iter()
        .enumerate()
        .map(|(i, &c)| c * offset.pow(&[i as u64]))
        .collect::<Vec<_>>();
    let composite_polynomial_evaluations: Vec<Fr> = extended_domain.fft(&shifted_coeffs);
    let (cp_merkle_tree, cp_commitment) = compute_merkle_root(&composite_polynomial_evaluations);
    (
        composite_polynomial_evaluations,
        cp_merkle_tree,
        cp_commitment,
    )
}

fn fri_layer(evals: &[Fr], beta: Fr, domain: &[Fr]) -> (Vec<Fr>, Vec<Fr>) {
    let n = evals.len();
    let half = n / 2;

    let mut next_evals = Vec::with_capacity(half);
    let mut next_domain = Vec::with_capacity(half);

    for i in 0..half {
        let x = domain[i];
        let fx = evals[i];
        let fnx = evals[i + half];

        let two_inv = Fr::from(2u64).inverse().unwrap();
        let g = (fx + fnx) * two_inv;
        let h = (fx - fnx) * two_inv * x.inverse().unwrap();

        next_evals.push(g + beta * h);
        next_domain.push(x * x);
    }
    (next_evals, next_domain)
}

fn calculate_fri_commitments(
    extended_evals: &[Fr],
    extended_domain: &[Fr],
    extended_trace_commitment: &[u8; 32],
    cp_commitment: &[u8; 32],
) -> (Vec<Vec<Vec<[u8; 32]>>>, Vec<[u8; 32]>, Vec<Vec<Fr>>) {
    let mut fiat_shamir_seed: Vec<u8> = Vec::new();
    fiat_shamir_seed.extend_from_slice(extended_trace_commitment);
    fiat_shamir_seed.extend_from_slice(cp_commitment);
    let mut current_evals = extended_evals.to_vec();
    let mut current_domain = extended_domain.to_vec();
    let mut commitments: Vec<[u8; 32]> = Vec::new();
    let mut merkle_trees: Vec<Vec<Vec<[u8; 32]>>> = Vec::new();
    let mut evals: Vec<Vec<Fr>> = Vec::new();
    while current_evals.len() > 1 {
        // add the latest committment to the seed
        if !commitments.is_empty() {
            let last_commitment = commitments.last().cloned().unwrap();
            fiat_shamir_seed.extend_from_slice(&last_commitment);
        }
        let beta = random_fr_from_hash(&fiat_shamir_seed);
        let (next_evals, next_domain) = fri_layer(&current_evals, beta, &current_domain);
        let (next_merkle_tree, next_commitment) = compute_merkle_root(&next_evals);
        commitments.push(next_commitment);
        merkle_trees.push(next_merkle_tree);
        evals.push(next_evals.clone());
        current_evals = next_evals;
        current_domain = next_domain;
    }
    (merkle_trees, commitments, evals)
}

fn extract_path_from_merkle_tree(index: &usize, merkle_tree: &[Vec<[u8; 32]>]) -> Vec<[u8; 32]> {
    let mut path = Vec::with_capacity(merkle_tree.len());
    let mut current_index = index.clone();
    for level in merkle_tree.iter().take(merkle_tree.len().saturating_sub(1)) {
        // push the sibling of the current node, since the verifier will need the siblings to recalculate the root
        let sibling_index = if current_index % 2 == 0 {
            current_index + 1
        } else {
            current_index - 1
        };
        path.push(level[sibling_index]);

        current_index /= 2;
    }
    path
}

fn calculate_query_proof(
    q: usize,
    blowup_factor: usize,
    extended_trace_evals: &[Fr],
    extended_trace_merkle_tree: &[Vec<[u8; 32]>],
    cp_evals: &[Fr],
    cp_merkle_tree: &[Vec<[u8; 32]>],
    fri_evals: &[Vec<Fr>],
    fri_merkle_tree: &[Vec<Vec<[u8; 32]>>],
) -> QueryProof {
    let f_x = ValueAndPath {
        value: extended_trace_evals[q],
        path: extract_path_from_merkle_tree(&q, extended_trace_merkle_tree),
    };
    let gq: usize = (q + blowup_factor) % extended_trace_evals.len();
    let f_gx = ValueAndPath {
        value: extended_trace_evals[gq],
        path: extract_path_from_merkle_tree(&gq, extended_trace_merkle_tree),
    };
    let ggq = (gq + blowup_factor) % extended_trace_evals.len();
    let f_g2x = ValueAndPath {
        value: extended_trace_evals[ggq],
        path: extract_path_from_merkle_tree(&ggq, extended_trace_merkle_tree),
    };
    let mut cp_pairs = Vec::with_capacity(1 + fri_evals.len());
    cp_pairs.push(CpPair {
        cp_x: ValueAndPath {
            value: cp_evals[q],
            path: extract_path_from_merkle_tree(&q, cp_merkle_tree),
        },
        cp_minus_x: ValueAndPath {
            value: cp_evals[(cp_evals.len() / 2 + q) % cp_evals.len()],
            path: extract_path_from_merkle_tree(
                &((cp_evals.len() / 2 + q) % cp_evals.len()),
                cp_merkle_tree,
            ),
        },
    });

    for (sub_cp_evals, sub_cp_merkle_tree) in fri_evals.iter().zip(fri_merkle_tree.iter()) {
        let index = q % sub_cp_evals.len();
        let negative_index = (sub_cp_evals.len() / 2 + index) % sub_cp_evals.len();
        cp_pairs.push(CpPair {
            cp_x: ValueAndPath {
                value: sub_cp_evals[index],
                path: extract_path_from_merkle_tree(&index, sub_cp_merkle_tree),
            },
            cp_minus_x: ValueAndPath {
                value: sub_cp_evals[negative_index],
                path: extract_path_from_merkle_tree(&negative_index, sub_cp_merkle_tree),
            },
        });
    }

    QueryProof {
        f_x,
        f_gx,
        f_g2x,
        cp_pairs,
    }
}

fn calculate_query_proofs(
    num_queries: usize,
    blowup_factor: usize,
    extended_trace_evals: &[Fr],
    extended_trace_merkle_tree: &[Vec<[u8; 32]>],
    cp_evals: &[Fr],
    cp_merkle_tree: &[Vec<[u8; 32]>],
    fri_evals: &[Vec<Fr>],
    fri_merkle_tree: &[Vec<Vec<[u8; 32]>>],
) -> Vec<QueryProof> {
    let mut fiat_shamir_seed: Vec<u8> = Vec::new();
    fiat_shamir_seed.extend_from_slice(extended_trace_merkle_tree.last().unwrap().last().unwrap());
    fiat_shamir_seed.extend_from_slice(cp_merkle_tree.last().unwrap().last().unwrap());
    for commitment in fri_merkle_tree.iter().flat_map(|tree| tree.last()) {
        for node in commitment {
            fiat_shamir_seed.extend_from_slice(node);
        }
    }

    let mut query_proofs = Vec::with_capacity(num_queries);

    for _ in 0..num_queries {
        let q = integer_from_hash(&fiat_shamir_seed, extended_trace_evals.len());
        let query_proof = calculate_query_proof(
            q,
            blowup_factor,
            extended_trace_evals,
            extended_trace_merkle_tree,
            cp_evals,
            cp_merkle_tree,
            fri_evals,
            fri_merkle_tree,
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
        query_proofs.push(query_proof);
    }

    query_proofs
}

pub fn generate_prover_values(
    a0: Fr,
    a1: Fr,
    target_index: usize,
    blowup_factor: usize,
    num_queries: usize,
) -> ProverValues {
    let trace_length = next_power_of_two(target_index + 1);
    let offset = Fr::GENERATOR;
    let trace = calculate_fib_square_trace(a0, a1, trace_length);
    let (trace_polynomial, domain) = calculate_trace_polynomial(&trace);
    let (extended_trace, extended_domain) =
        lde_on_coset(&trace_polynomial, trace_length, blowup_factor, offset);
    let (extended_trace_merke_tree, extended_trace_commitment) =
        compute_merkle_root(&extended_trace);
    let (cp_evaluations, cp_merkle_tree, cp_commitment) = calculate_composite_polynomial(
        &trace_polynomial,
        &domain,
        &extended_domain,
        offset,
        *trace.get(target_index).unwrap(),
        target_index,
        trace_length,
    );
    let extended_domain_points: Vec<Fr> = extended_domain.elements().map(|x| x * offset).collect();
    let (fri_merkle_trees, fri_commitments, fri_evals) = calculate_fri_commitments(
        &cp_evaluations,
        &extended_domain_points,
        &extended_trace_commitment,
        &cp_commitment,
    );
    let query_proofs = calculate_query_proofs(
        num_queries,
        blowup_factor,
        &extended_trace,
        &extended_trace_merke_tree,
        &cp_evaluations,
        &cp_merkle_tree,
        &fri_evals,
        &fri_merkle_trees,
    );
    ProverValues {
        extended_trace_commitment,
        composite_polynomial_commitment: cp_commitment,
        fri_commitments,
        queries_proofs: query_proofs,
    }
}
