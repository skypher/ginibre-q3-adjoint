# Exact AM-GM advance for the rank-six orbit ring `O_11`

Date: 2026-07-24

## Status

This note proves the first previously unresolved all-exponent chamber in the
rank-six odd simple-current orbit ring `O_11` and supplies exact certificates
for 1,583 further residual orthants.

The complete `O_11` theorem is not yet proved. The new results are:

1. the chamber

   ```text
   (B_1^-)^(2+2p) (B_5^+)^(1+2q),       p,q>=0,
   ```

   is nonnegative for all exponents;
2. 25 further residual keys have hand-audited exact AM-GM certificates;
3. a deterministic ledger supplies exact AM-GM certificates for another
   1,558 residual keys in 294 support/parity chambers.

By odd simple-current lifting, completion of the remaining `O_11` chambers
would prove full `GKS2*` and the whole partial-character column in
`SU(2)_11`.

## 1. Exact algebraic model

The six trace nodes of `O_11` are the six real roots of

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

Every verifier constructs a rational Sturm sequence, proves that six rational
intervals isolate all six roots, evaluates the orbit characters with rational
interval arithmetic, and accepts an inequality only when the interval
endpoints prove it. Floating-point allocation searches are not proof evidence.

## 2. Capacitated weighted AM-GM certificate

After fixing a support/sign/parity chamber and absorbing the residual-orthant
floor into the coefficients, write

```text
F(r)=sum_(j in P) C_j product_l lambda_(j,l)^r_l
     -sum_(n in N) D_n product_l mu_(n,l)^r_l.
```

For each negative term choose rational weights `alpha_(n,j)>=0` satisfying

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l),
sum_n alpha_(n,j)D_n <= C_j.
```

Weighted AM-GM then gives `F(r)>=0` for every nonnegative integral residual
vector. All certificates here use integer weights with denominator 100. The
verifiers check products of algebraic intervals raised to integer powers; no
logarithm or approximate convex comparison is trusted.

## 3. Complete first-chamber theorem

`verify_su2_o11_first_chamber_exact.cpp` partitions the lattice into

```text
q>=4, p>=0;
q=3,  p>=0;
q=2,  p>=1;
q=1,  p>=1;
q=0,  p>=2.
```

Each region has an exact denominator-100 AM-GM certificate. The four omitted
points

```text
(0,0), (1,0), (0,1), (0,2)
```

are evaluated by integral orbit-fusion dynamic programming and all have
corner zero. Therefore

```text
[V_0 tensor V_0]
 (B_1 tensor 1-1 tensor B_1)^(2+2p)
 (B_5 tensor 1+1 tensor B_5)^(1+2q) >= 0
```

for every `p,q>=0`.

The strict replay reports

```text
SU2_O11_FIRST_CHAMBER_EXACT PASS
roots=6 spectral_pairs=15 regions=5
amgm_denominator=100 exact_leaves=4 threads=5
```

This replaces the earlier bounded `101 x 101` scan by an all-exponent proof.

## 4. Further exact regions

`verify_su2_o11_amgm_frontier_exact.cpp` proves 25 residual keys with 26
certificates; one key is divided into a boundary strip and a translated tail.
Its strict replay is

```text
SU2_O11_AMGM_FRONTIER_EXACT PASS
certificates=26 residual_keys=25 denominator=100 threads=4
```

A four-worker mixed-integer search then generated denominator-100 candidate
allocations for a larger conservative census. The numerical search is not
trusted. The exact selector payload

```text
certificates/su2_o11_amgm_batch_ledger.b64
```

contains only support, parity, residual code, spectral-pair selectors, and
integer weights. Its SHA-256 digest is

```text
4bad077b05a9f5554729da30af081728abe1578b70dc271a50813dab01794c69.
```

`verify_su2_o11_amgm_batch_exact.cpp` independently reconstructs every
chamber and checks:

1. every negative spectral pair occurs exactly once;
2. every allocation row sums to its denominator;
3. every geometric-base inequality holds in each residual coordinate;
4. every positive-capacity inequality holds;
5. all `(support,parity,residual)` keys are distinct.

The decoded ledger has 1,558 records in 294 support/parity chambers. Exact
text-ledger replay passed all records in sixteen ranges:

```text
15 ranges: 100 records PASS
 1 range :  58 records PASS
total   : 1558 records PASS.
```

The compact binary decoder and representative range replay were checked
independently.

## 5. Remaining frontier

These results do not prove the complete rank-six orbit theorem. A conservative
reduced-budget transport census left 142 candidate regions for which one
global denominator-100 AM-GM allocation exceeds positive capacity. This is
not evidence of negativity. Those regions require lattice splitting,
translation with exact face certificates, a different transport, or a
structural Turan argument. An authoritative full 100-ray/100-shift census is
also still needed.

## 6. Reproduction

Build the compact exact verifiers with strict warnings:

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

Build the larger ledger verifier with `-O2`. The source-matched transcript is
`certificates/su2_o11_amgm_exact.log`.
