use ark_bn254::Fr;
use ark_ff::{FftField, Field};
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DensePolynomial,
};

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
    trace: Vec<Fr>,
) -> (DensePolynomial<Fr>, GeneralEvaluationDomain<Fr>) {
    let domain = GeneralEvaluationDomain::<Fr>::new(trace.len()).unwrap();

    let evaluations = Evaluations::from_vec_and_domain(trace, domain);

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

pub fn generate_prover_values(a0: Fr, a1: Fr, trace_length: usize) -> ProverValues {
    let trace = calculate_fib_square_trace(a0, a1, trace_length);
    let (trace_polynomial, _domain) = calculate_trace_polynomial(trace);
    let (extended_trace, _) = lde_on_coset(&trace_polynomial, trace_length, 8, Fr::GENERATOR);

    ProverValues { extended_trace }
}
