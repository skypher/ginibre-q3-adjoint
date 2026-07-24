# Exact AM-GM advance for the rank-six orbit ring `O_11`

Date: 2026-07-24

## Status

This note proves the first previously unresolved all-exponent chamber in the
rank-six odd simple-current orbit ring `O_11` and supplies exact certificates
for 1,587 further residual orthants.

The complete `O_11` theorem is not yet proved.  The new results are:

1. the two-label chamber

   ```text
   (B_1^-)^(2+2p) (B_5^+)^(1+2q),       p,q>=0,
   ```

   is nonnegative for all exponents;
2. 29 further residual keys have hand-audited exact AM-GM certificates;
3. a deterministic ledger supplies exact AM-GM certificates for another
   1,558 residual keys in 294 support/parity chambers.

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
only when the rational interval endpoints prove it. Floating-point output
from the allocation search is not part of the proof.

## 2. Capacitated weighted AM-GM lemma

Fix a support/sign/parity chamber and a residual orthant. After absorbing the
orthant floor into the coefficients, its corner has the form

```text
F(r)=sum_(j in P) C_j product_l lambda_(j,l)^r_l
     -sum_(n in N) D_n product_l mu_(n,l)^r_l,
```

where every coefficient and base is positive.

**Lemma 2.1.** Suppose that for each negative term `n` there are rational
weights `alpha_(n,j)>=0` such that

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l)
                                      for every residual coordinate l,
sum_n alpha_(n,j)D_n <= C_j           for every positive term j.
```

Then `F(r)>=0` for every nonnegative integral residual vector `r`.

**Proof.** The base inequalities and weighted AM-GM give

```text
D_n product_l mu_(n,l)^r_l
 <= sum_j alpha_(n,j)D_n product_l lambda_(j,l)^r_l.
```

Sum over `n` and use the capacity inequalities. QED.

All certificates below use integer weights with denominator one hundred. The
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

Each region has a denominator-100 AM-GM certificate. The four omitted points

```text
(p,q)=(0,0),(1,0),(0,1),(0,2)
```

are evaluated by exact integral orbit-fusion dynamic programming and all have
corner zero. Therefore

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

## 4. Hand-audited frontier block

`verify_su2_o11_amgm_frontier_exact.cpp` proves 29 distinct residual keys by
30 certificates; one key is divided into a boundary strip and a translated
tail. The block covers support codes `86,89,95,98,101` and related immediately
adjacent regimes.

Its strict replay is

```text
SU2_O11_AMGM_FRONTIER_EXACT PASS
certificates=30 residual_keys=29 denominator=100 threads=4
```

## 5. Batch exact ledger

A reduced-budget parallel transport census was used only to generate candidate
regions. After the complete first chamber and the hand frontier block were
removed, it flagged 1,700 residual keys not covered by that diagnostic's
coordinatewise Hall, braid, rational-fan, or shifted-face routines.

A four-worker mixed-integer search found denominator-100 AM-GM allocations
with capacity at most one for 1,558 of those keys. The numerical search is
not trusted. The file

```text
certificates/su2_o11_amgm_batch_ledger.b64
```

contains only the discrete selector data: support, parity, residual code, and
integer allocation weights. Its SHA-256 digest is

```text
4bad077b05a9f5554729da30af081728abe1578b70dc271a50813dab01794c69.
```

`verify_su2_o11_amgm_batch_exact.cpp` independently reconstructs every
chamber and checks, with rational algebraic intervals:

1. every negative spectral pair occurs exactly once;
2. each allocation row sums to its stated denominator;
3. every geometric base inequality holds in every active residual coordinate;
4. every positive capacity inequality holds;
5. all `(support,parity,residual)` keys are distinct.

The decoded ledger has 1,558 records and spans 294 support/parity chambers.
For bounded wall-time and better cache locality, the recorded replay selects
sixteen consecutive ranges of at most one hundred records. The exact C++
verifier audits the full key set before checking a selected range. The range
replays collectively report

```text
15 ranges: certificates=100 residual_keys=100 PASS
 1 range : certificates=58  residual_keys=58  PASS
 total   : certificates=1558 residual_keys=1558 PASS.
```

The verifier supports explicit range selection and one to five worker threads.
The text-ledger replay completed all 1,558 exact records; the compact binary
decoder was independently checked on representative ranges.

## 6. Remaining frontier

The 1,558 certificates settle every region in the batch ledger, but they do
not prove the complete rank-six orbit theorem. There remain 142 regions from
the conservative census for which a single global denominator-100 AM-GM
allocation exceeds positive capacity. This does not indicate negativity: it
only means those regions need a lattice split, translation with exact face
certificates, a different transport, or a structural Turan argument.

The next useful task is to classify the 142 residual regions by lattice faces
and seek endpoint-lattice or translated-tail certificates, while separately
completing an authoritative full 100-ray/100-shift census.

## 7. Reproduction

Build the two compact exact verifiers with

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

The larger ledger verifier is faster to compile with `-O2`:

```text
g++ -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
  -Wsign-conversion -Wshadow -Werror -pthread \
  verify_su2_o11_amgm_batch_exact.cpp \
  -o verify_su2_o11_amgm_batch_exact
```

The source-matched transcript is
`certificates/su2_o11_amgm_exact.log`.
