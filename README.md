# STARK 101 - Rust Implementation

This repo contains my Rust implementation of a STARK prover from scratch as showcased in the S[TARK 101 hands-on tutorial](https://starkware.co/stark-101/)

This STARK will be built part by part as shown in the tutorial.

The goal of this STARK is to prove that the prover knows a term $a_1$ such that $a_{15} = 20058280215495444632052566758236617048289674862308296983290231865868158747890$ and $a_{0} = 1$, where $a_{n+2} = a_{n+1}^2 + a_n^2$.

##  Part1

In this part, the code contains:

- the Trace generation: $a_0, a_1, a_2, ..., a_{15}$.
- its polynomial interpolation: we find a polynomial $P$ such that $P(g^i) = a_i$ for $i=0,1,...,15$, where $g$ is a generator of the 2-aditive subgroup of size 16 of the field.
- and then its Low Degree Extension: we evaluate the polynomial $P$ on a larger set of points, which is the 2-aditive subgroup of size 64.
