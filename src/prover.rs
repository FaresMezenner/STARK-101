use ark_bn254::Fr;
use ark_ff::{FftField, Field, PrimeField, UniformRand, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DenseOrSparsePolynomial,
    univariate::DensePolynomial,
};
use sha2::Digest;

use crate::core::compute_merkle_root;
use crate::data::ProverValues;

fn calculate_fib_square_trace(a0: Fr, a1: Fr, length: usize) -> Vec<Fr> {
    let mut trace: Vec<Fr> = vec![a0, a1];

    for i in 2..length {
        let next = trace[i - 1].square() + trace[i - 2].square();
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
    target: Fr,
) -> (Vec<Fr>, [u8; 32]) {
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

    // calculating p1 = (f(X) - a_{15}) / (X - g^15)
    let p1_numerator: DenseOrSparsePolynomial<Fr> =
        DenseOrSparsePolynomial::from(f - &DensePolynomial::from_coefficients_vec(vec![target]));
    let p1_denominator: DenseOrSparsePolynomial<Fr> =
        DenseOrSparsePolynomial::from(DensePolynomial::from_coefficients_vec(vec![
            -g.pow(&[15u64]),
            Fr::from(1u64),
        ]));
    let (p1, p1_remainder) = p1_numerator.divide_with_q_and_r(&p1_denominator).unwrap();
    assert!(p1_remainder.is_zero(), "p1 division must be exact");

    // calculating p3 = (f((g^2)X) - f(gX)^2 - f(X)^2) * (X - g^15) * (X - g^14) / (X^16 - 1)
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

    let x_minus_g15 =
        DensePolynomial::from_coefficients_vec(vec![-g.pow(&[15u64]), Fr::from(1u64)]);
    let x_minus_g14 =
        DensePolynomial::from_coefficients_vec(vec![-g.pow(&[14u64]), Fr::from(1u64)]);
    let p3_numerator =
        DenseOrSparsePolynomial::from(&(&transition_numerator * &x_minus_g15) * &x_minus_g14);

    let p3_denominator = DenseOrSparsePolynomial::from(DensePolynomial::from_coefficients_vec(
        (0..=16)
            .map(|i| {
                if i == 0 {
                    -Fr::from(1u64)
                } else if i == 16 {
                    Fr::from(1u64)
                } else {
                    Fr::zero()
                }
            })
            .collect(),
    ));
    let (p3, p3_remainder) = p3_numerator.divide_with_q_and_r(&p3_denominator).unwrap();
    assert!(p3_remainder.is_zero(), "p3 division must be exact");

    // now calculating CP using p0, p1, p3 and random linear combination
    let mut rng = rand::thread_rng();
    let alpha_0 = Fr::rand(&mut rng);
    let alpha_1 = Fr::rand(&mut rng);
    let alpha_2 = Fr::rand(&mut rng);
    let composite_polynomial = &(&(&p0 * alpha_0) + &(&p1 * alpha_1)) + &(&p3 * alpha_2);
    let composite_polynomial_evaluations = extended_domain.fft(&composite_polynomial.coeffs);
    let cp_commitment = compute_merkle_root(&composite_polynomial_evaluations);
    (composite_polynomial_evaluations, cp_commitment)
}

fn random_fr_from_hash(input: &[u8]) -> Fr {
    let hash = sha2::Sha256::digest(input);
    let mut bytes = [0u8; 32];
    bytes.copy_from_slice(&hash);
    Fr::from_le_bytes_mod_order(&bytes)
}

fn fri_layer(evals: &[Fr], beta: Fr, domain: &[Fr]) -> (Vec<Fr>, Vec<Fr>) {
    let n = evals.len();
    let half = n / 2;

    let mut next_evals = Vec::with_capacity(half);
    let mut next_domain = Vec::with_capacity(half);

    for i in 0..half {
        let x = domain[i];
        let neg_x = domain[i + half];
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
) -> Vec<[u8; 32]> {
    let mut fiat_shamir_seed: Vec<u8> = Vec::with_capacity(64);
    fiat_shamir_seed.extend_from_slice(extended_trace_commitment);
    fiat_shamir_seed.extend_from_slice(cp_commitment);
    let mut current_evals = extended_evals.to_vec();
    let mut current_domain = extended_domain.to_vec();
    let mut commitments: Vec<[u8; 32]> = Vec::new();
    while current_evals.len() > 1 {
        // add the latest committment to the seed
        let last_commitment = commitments.last().cloned().unwrap_or([0u8; 32]);
        fiat_shamir_seed.extend_from_slice(&last_commitment);
        let beta = random_fr_from_hash(&fiat_shamir_seed);
        let (next_evals, next_domain) = fri_layer(&current_evals, beta, &current_domain);
        let next_commitment = compute_merkle_root(&next_evals);
        commitments.push(next_commitment);
        current_evals = next_evals;
        current_domain = next_domain;
    }
    commitments
}

pub fn generate_prover_values(a0: Fr, a1: Fr, trace_length: usize) -> ProverValues {
    let trace = calculate_fib_square_trace(a0, a1, trace_length);
    let (trace_polynomial, domain) = calculate_trace_polynomial(&trace);
    let (extended_trace, extended_domain) =
        lde_on_coset(&trace_polynomial, trace_length, 8, Fr::GENERATOR);
    let extended_trace_commitment = compute_merkle_root(&extended_trace);
    let (cp_evaluations, cp_commitment) = calculate_composite_polynomial(
        &trace_polynomial,
        &domain,
        &extended_domain,
        *trace.last().unwrap(),
    );
    let extended_domain_points: Vec<Fr> = extended_domain.elements().collect();
    let fri_commitments = calculate_fri_commitments(
        &cp_evaluations,
        &extended_domain_points,
        &extended_trace_commitment,
        &cp_commitment,
    );
    ProverValues {
        extended_trace_commitment,
        composite_polynomial_commitment: cp_commitment,
        fri_commitments,
    }
}
