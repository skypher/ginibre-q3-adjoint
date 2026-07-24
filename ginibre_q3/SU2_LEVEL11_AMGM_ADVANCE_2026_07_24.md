# Exact AM-GM advance for the rank-six orbit ring `O_11`

Date: 2026-07-24

## Status

This note proves the first previously unresolved all-exponent chamber in the
rank-six odd simple-current orbit ring `O_11` and supplies exact certificates
for 25 further residual orthants.

The complete `O_11` theorem is not yet proved.  The new results are:

1. the two-label chamber

   ```text
   (B_1^-)^(2+2p) (B_5^+)^(1+2q),       p,q>=0,
   ```

   is nonnegative for all exponents;
2. 25 further residual keys have hand-audited exact AM-GM certificates,
   using 26 certificates because one orthant is split into a strip and tail.

By odd simple-current lifting, completion of the remaining `O_11` chambers
would prove full `GKS2*` and the whole partial-character column in
`SU(2)_11`.

## 1. Spectral arithmetic

The six trace nodes of `O_11` are represented by the six real roots of

```text
f(x)=x^6-5x^5+5x^4+6x^3-7x^2-2x+1.
```

In the even-lift orbit basis,

```text
B_1=x,
B_2=x^2-x-1,
B_3=x^3-2x^2-x+1,
B_4=x^4-3x^3+3x,
B_5=x^5-4x^4+2x^3+5x^2-2x-1,
```

and the trace weight is `(3-x)/13`.

Every exact verifier in this note constructs a rational Sturm sequence,
checks that six rational intervals isolate all six roots, evaluates all orbit
characters by rational interval arithmetic, and accepts a sign or inequality
only when the rational interval endpoints prove it.  Floating-point output
from the allocation search is not part of the proof.

## 2. Capacitated weighted AM-GM lemma

Fix a support/sign/parity chamber and a residual orthant.  After absorbing the
orthant floor into the coefficients, its corner has the form

```text
F(r)=sum_(j in P) C_j product_l lambda_(j,l)^r_l
     -sum_(n in N) D_n product_l mu_(n,l)^r_l,
```

where every coefficient and base is positive.

**Lemma 2.1.**  Suppose that for each negative term `n` there are rational
weights `alpha_(n,j)>=0` such that

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l)
                                      for every residual coordinate l,
sum_n alpha_(n,j)D_n <= C_j           for every positive term j.
```

Then `F(r)>=0` for every nonnegative integral residual vector `r`.

**Proof.**  The base inequalities and weighted AM-GM give

```text
D_n product_l mu_(n,l)^r_l
 <= sum_j alpha_(n,j)D_n product_l lambda_(j,l)^r_l.
```

Sum over `n` and use the capacity inequalities.  QED.

All certificates below use integer weights with denominator one hundred.  The
verifiers raise the algebraic base intervals to the corresponding integer
powers, so no logarithm or approximate convex comparison is trusted.

## 3. Complete proof of the first chamber

For

```text
W_(p,q)=(B_1 tensor 1-1 tensor B_1)^(2+2p)
        (B_5 tensor 1+1 tensor B_5)^(1+2q),
```

`verify_su2_o11_first_chamber_exact.cpp` partitions the nonnegative lattice
into five regions:

```text
q>=4,                         p>=0,
q=3,                          p>=0,
q=2,                          p>=1,
q=1,                          p>=1,
q=0,                          p>=2.
```

Each region has a denominator-100 AM-GM certificate.  The four omitted points

```text
(p,q)=(0,0),(1,0),(0,1),(0,2)
```

are evaluated by exact integral orbit-fusion dynamic programming and all have
corner zero.  Therefore

```text
[V_0 tensor V_0] W_(p,q) >= 0       for every p,q>=0.
```

The strict replay reports

```text
SU2_O11_FIRST_CHAMBER_EXACT PASS
roots=6 spectral_pairs=15 regions=5
amgm_denominator=100 exact_leaves=4 threads=5
```

This replaces the earlier bounded `101 x 101` scan by an all-exponent proof.

## 4. Exact frontier block

`verify_su2_o11_amgm_frontier_exact.cpp` proves 25 distinct residual keys by
26 certificates; one key is divided into a boundary strip and a translated
tail. The strict replay reports

```text
SU2_O11_AMGM_FRONTIER_EXACT PASS
certificates=26 residual_keys=25 denominator=100 threads=4
```

The block includes two-, three-, and four-variable residual interiors and
shows that the weighted AM-GM method is not restricted to the original
two-label chamber. Every certificate is reconstructed from the six algebraic
trace nodes and checked with rational interval arithmetic.

## 5. Remaining frontier

The complete rank-six orbit theorem is still open. A reduced-budget transport
census and parallel mixed-integer allocation search produced a much larger
candidate ledger, but that exploratory batch is not included in this theorem:
its compact repository replay still needs a complete source-matched run.

The next rigorous task is to extend the exact frontier ledger with translated
orthants and exact boundary faces, and then run an authoritative full
100-ray/100-shift census.

## 6. Reproduction

Build the two exact verifiers with

```text
g++ -O3 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
  -Wsign-conversion -Wshadow -Werror -pthread \
  verify_su2_o11_first_chamber_exact.cpp \
  -o verify_su2_o11_first_chamber_exact

g++ -O3 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
  -Wsign-conversion -Wshadow -Werror -pthread \
  verify_su2_o11_amgm_frontier_exact.cpp \
  -o verify_su2_o11_amgm_frontier_exact
```


The full source-matched transcript is
`certificates/su2_o11_amgm_exact.log`.
