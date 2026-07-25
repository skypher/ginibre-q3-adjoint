# Exact residual-dimension `<=2` `GKS2*` for the rank-seven orbit ring `O_13`

Date: 2026-07-25

## Result

Let `O_13` be the rank-seven simple-current orbit ring of `SU(2)_13`, with
even-lift basis

```text
B_0=V_0, B_1=V_2, ..., B_6=V_12.
```

After support-overlap reduction, fix the sign and parity of every active label
and write each active exponent as

```text
p_a=p_a^0+2r_a,      r_a>=0,
```

where `p_a^0` is one or two.  The **residual dimension** is the number of
coordinates `r_a` allowed to vary.

**Theorem.** Every support-disjoint even-minus scalar chamber in `O_13` of
residual dimension at most two has nonnegative scalar corner for all residual
exponents.  By reversing support-overlap reduction, the same holds without
the support-disjoint assumption.  Adjoining one plus factor therefore also
gives the corresponding partial-character coefficients.

This is an all-exponent theorem, not a bounded scan.

## 1. Exact low-dimensional census

The exact rational-Sturm/direct-Hall census gives:

```text
residual dimension 0
  total regimes                   7,448
  pointwise positive                665
  exact direct Hall               6,669
  exact integral leaves             114
  unresolved                           0

residual dimension 1
  total regimes                  36,042
  pointwise positive              2,724
  exact direct Hall              32,808
  post-Hall regimes                 510

residual dimension 2
  total regimes                  72,570
  pointwise positive              4,620
  exact direct Hall              63,914
  post-Hall regimes               4,036
```

Thus the theorem covers

```text
7,448+36,042+72,570=116,060
```

residual regimes.

The census isolates all seven trace nodes by rational Sturm sequences,
classifies exact cyclotomic zero factors in `Q(zeta_30)`, and uses outward
MPFR/logarithmic intervals only after the algebraic signs and zeros have been
settled exactly.

## 2. One-variable post-Hall sector

Of the 510 exact post-Hall rays:

```text
unique dominant spectral pair + finite exact leaves     508
separated-pole theorem, Lemma 22H3S                        2
```

For each of the 508 dominant-pair rays, one positive spectral pair has a
strictly larger squared residual base than every negative pair.  Once its
mass dominates the total negative mass it remains dominant forever.  Exact
integral fusion arithmetic checks every earlier exponent.  The largest tail
threshold is 23 and the total number of exact leaves is 3,142.

The two additional rays are the coordinate faces of

```text
(B_1^-)^(1+2p) (B_6^-)^(1+2q),
```

which is the separated-pole/top-radius chamber already proved for all indices
by Lemma 22H3S.

## 3. Two-variable post-Hall sector

The exact partition of all 4,036 post-Hall keys is:

```text
single dominant pair + finite simplex                 411
rational weighted AM-GM                               255
exact ordered-braid/Farey transport                 2,780
uniform-margin replacement for shifted transport      272
uniform-margin batch for final hard set                314
two individually certified margin chambers              2
top-orbit identity and tadpole-ray theorem                1
separated-pole theorem                                    1
                                                     -----
total                                                4,036
```

The key sets are pairwise disjoint and their literal union equals the exact
post-Hall failure set.

### Weighted AM-GM certificates

For a negative spectral monomial `D_n mu_n^r`, rational weights satisfy

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l),
sum_n alpha_(n,j)D_n <= C_j.
```

The replay raises algebraic intervals to integer powers or compares their
logs with directed rounding; numerical optimization is used only to propose
integer weights and is never trusted by the verifier.

### Uniform-margin certificates

For the translated and final hard sectors, the convex allocation has a
strict logarithmic advantage `delta>0` in both residual coordinates.  A
finite capacity-overbooking factor is therefore killed once `p+q` exceeds a
certified threshold.  Exact integral fusion arithmetic checks the remaining
finite simplex.

The 272 shifted cases require 328,711 exact finite values; the 314 final
batch requires 548,849.  The two individually certified chambers require 120
and 1,770 values.

### Top-orbit reduction

In `O_13`,

```text
B_6^2=B_0+B_1,
D_1=D_6 S_6,
```

where `D_a=B_a tensor 1-1 tensor B_a` and
`S_a=B_a tensor 1+1 tensor B_a`.  Hence key `(245,2,3)` reduces exactly to a
single top-orbit/tadpole ray, covered by Lemma 22H3S.  The replay includes an
integral matrix check of both identities.

## 4. Reproducibility

The deterministic replay archive is stored as the text chunks

```text
certificates/su2_o13_lowdim_exact_replay/su2_o13_lowdim_exact_replay.tar.gz.b64.part*
```

with SHA-256

```text
7497005a39233054e35efe3449ae8df8c7693f0f9def96a2b3b18eb7c8cd5796
```

Reconstruct, extract, and run it with

```sh
cd certificates/su2_o13_lowdim_exact_replay
cat su2_o13_lowdim_exact_replay.tar.gz.b64.part* \
  | base64 -d > su2_o13_lowdim_exact_replay.tar.gz
sha256sum -c SHA256SUMS
tar -xzf su2_o13_lowdim_exact_replay.tar.gz
./replay_o13_lowdim/replay.sh
```

on a system with C++20, Boost.Multiprecision, GMP, the MPFR 4 runtime ABI, and
Python 3.  The replay uses five worker threads for the expensive C++ checks.

The clean extracted replay reports, among other lines,

```text
O13_LOWDIM_DIRECT_FAST_EXACT PASS chambers=7448 threads=5
O13_DIM1_TAIL_EXACT PASS keys=508 maximum_threshold=23 exact_leaves=3142
O13_DIM2_TAIL_EXACT PASS keys=411 maximum_threshold=20 exact_simplex_leaves=7728
O13_DIM2_AMGM_EXACT PASS certificates=255
O13_DIM2_BRAID_FAN_FAST_EXACT PASS braid=1480 fan=0
O13_DIM2_BRAID_FAN_FAST_EXACT PASS braid=0 fan=1300
O13_DIM2_MARGIN_GENERIC_EXACT PASS certificates=314 finite_checks=548849
O13_DIM2_MARGIN_GENERIC_EXACT PASS certificates=272 finite_checks=328711
O13_LOWDIM_KEY_AUDIT PASS dim1_post_hall=510 dim2_post_hall=4036
O13_LOWDIM_REPLAY PASS
```

## 5. Remaining frontier

Full scalar `GKS2*` for `O_13`, and hence full `SU(2)_13`, remains open.
The next unresolved sector has residual dimension three.  The diagnostic
rank-seven census contains 11,486 post-Hall three-variable regimes before
applying higher-dimensional braid, fan, AM-GM, symmetry, and shifted-face
certificates.
