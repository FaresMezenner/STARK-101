# STARK 101 - Rust Implementation

This repo contains my Rust implementation of a STARK prover from scratch as showcased in the S[TARK 101 hands-on tutorial](https://starkware.co/stark-101/)

This STARK will be built part by part as shown in the tutorial.

The goal of this STARK is to prove that the prover knows a term $a_1$ such that $a_{15} = 20058280215495444632052566758236617048289674862308296983290231865868158747890$ and $a_{0} = 1$, where $a_{n+2} = a_{n+1}^2 + a_n^2$.

##  Part1

In this part, the code contains:

- the Trace generation: $a_0, a_1, a_2, ..., a_{15}$.
- its polynomial interpolation: we find a polynomial $f$ such that $f(g^i) = a_i$ for $i=0,1,...,15$, where $g$ is a generator of the 2-aditive subgroup of size 16 of the field.
- and then its Low Degree Extension: we evaluate the polynomial $f$ on a larger set of points, which is the 2-aditive subgroup of size 64.

##  Part 2

Now we will turn the statements to statements about polynomials as follows:

- $a_{0} = 1 \iff f(g^0) = 1 \iff f(g^0) - 1 = 0 \iff p_{0}(X) = \frac{f(X) - 1}{X - g^0} \text{ is a polynomial}$
- $a_{15} = 20058280215495444632052566758236617048289674862308296983290231865868158747890 \iff f(g^{15}) = 20058280215495444632052566758236617048289674862308296983290231865868158747890 \iff f(g^{15}) - 20058280215495444632052566758236617048289674862308296983290231865868158747890 = 0 \iff p_{2}(X) = \frac{f(X) - 20058280215495444632052566758236617048289674862308296983290231865868158747890}{X - g^{15}} \text{ is a polynomial}$
- $a_{n+2} = a_{n+1}^2 + a_{n}^2 \iff f(g^{n+2}) = f(g^{n+1})^2 + f(g^{n})^2 \iff f(g^{n+2}) - f(g^{n+1})^2 - f(g^{n})^2 = 0 \iff p_{3}(X) = \frac{f(X) - f(gX)^2 - f(g^2X)^2}{\frac{X^{16} -1}{(X-g^{15})(X-g^{14})}} \text{ is a polynomial}$

We will then calculate the Composition Polynomial: $CP = \alpha_0 p_0 + \alpha_1 p_1 + \alpha_2 p_2 + \alpha_3 p_3$, where $\alpha_i$ are random field elements. The prover will then commit to the coefficients of $CP$ using a Merkle tree.
Proving that CP is a polynomial implies that $p_0, p_1, p_2, p_3$ are polynomials, which in turn implies that the original statements about the trace are true.

##  Part 3

Instead of proving that $CP$ is a polynomial, we will prove that $CP$ is CLOSE to a polynomial of degree $d$ using FRI as follows:

- The prover has already committed to the evaluations of $CP$ on the extended domain (the 2-aditive subgroup of size 64).
- The verifier will then sample a random point $\beta$ and ask the prover to provide commitment to $CP'$ where $CP'(X) = g(X) + \beta·h(X)$ such that $CP(X) = g(X) + X·h(X)$, in other words, $g(X)$ is the even part of $CP$ and $h(X)$ is the odd part of $CP$.
- The commitment happens using the evaluation of $CP'$ on smaller and fresh domain.
- repeat the above steps until the bounded degree becomes $d<1$.
