use ark_bn254::Fr;

use crate::{prover::generate_prover_values, verifier::verify};

mod core;
mod data;
mod prover;
mod verifier;

fn main() {
    // prover calculates the extended trace polynomial
    // with these values, the last value in the trace is: 20058280215495444632052566758236617048289674862308296983290231865868158747890
    let a0 = Fr::from(1u64);
    let a1 = Fr::from(10u64);
    let a15: Fr = "20058280215495444632052566758236617048289674862308296983290231865868158747890"
        .parse()
        .expect("valid Fr element");
    let prover_values = generate_prover_values(a0, a1, 16, 8, 20);
    verify(prover_values, 16, 8, 20, a15);
}
