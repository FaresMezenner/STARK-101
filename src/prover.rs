use ark_bn254::Fr;
use ark_ff::{FftField, Field, UniformRand, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DenseOrSparsePolynomial,
    univariate::DensePolynomial,
};

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
) -> [u8; 32] {
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
    cp_commitment
}

pub fn generate_prover_values(a0: Fr, a1: Fr, trace_length: usize) -> ProverValues {
    let trace = calculate_fib_square_trace(a0, a1, trace_length);
    let (trace_polynomial, domain) = calculate_trace_polynomial(&trace);
    let (extended_trace, extended_domain) =
        lde_on_coset(&trace_polynomial, trace_length, 8, Fr::GENERATOR);
    let cp_commitment = calculate_composite_polynomial(
        &trace_polynomial,
        &domain,
        &extended_domain,
        *trace.last().unwrap(),
    );
    ProverValues {
        extended_trace_commitment: compute_merkle_root(&extended_trace),
        composite_polynomial_commitment: cp_commitment,
    }
}
