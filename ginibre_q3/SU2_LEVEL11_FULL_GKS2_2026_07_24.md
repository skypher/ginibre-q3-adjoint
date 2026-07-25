# Full `GKS2*` for `SU(2)_11`

Date: 2026-07-24

## Result

The rank-six simple-current orbit ring `O_11` satisfies scalar `GKS2*` for
every signed word of arbitrary length. By the odd-level simple-current lifting
theorem, the complete fusion ring `SU(2)_11` therefore satisfies scalar
`GKS2*`. Adjoining one plus factor proves that every partial-character
coefficient is nonnegative.

Equivalently, every even-minus signed word in the irreducible real character
basis of `SU(2)_11` has nonnegative vacuum corner, and for every target `c`,

```text
[V_c tensor V_0]
 product_i(V_(a_i) tensor V_0
           +epsilon_i V_0 tensor V_(a_i)) >= 0.
```

## 1. Exact spectral model

Use the even-lift orbit basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6, B_4=V_8, B_5=V_10.
```

The six trace nodes are the six real roots of

```text
x^6-5x^5+5x^4+6x^3-7x^2-2x+1.
```

The orbit characters are

```text
B_1=x,
B_2=x^2-x-1,
B_3=x^3-2x^2-x+1,
B_4=x^4-3x^3+3x,
B_5=x^5-4x^4+2x^3+5x^2-2x-1,
```

with trace weight `(3-x)/13`. Rational Sturm sequences certify one root in
each recorded decimal interval. All subsequent sign and capacity decisions use
outward interval arithmetic. Weighted geometric inequalities are replayed with
directed-rounding MPFR.

## 2. Exhaustive chamber census

Support overlap reduces to support-disjoint sign chambers. Fixing minimum
exponent parities and splitting each residual orthant by its zero coordinates
produces exactly

```text
27,962 residual regimes.
```

The exact exhaustive partition is

```text
pointwise nonnegative spectral sums       2,882
direct coordinatewise Hall transport    22,443
exact braid transport                       371
post-braid certificate archive            2,266
                                         ------
total                                    27,962
```

The fast universal census validates the rational root brackets first and then
uses outward-rounded interval operations. Its exact output is

```text
O11_FULL_CENSUS_FAST_EXACT PASS
all=27962 pointwise=2882 direct=22443 braid=371 post_braid=2266 threads=5
```

## 3. The 2,266 post-braid regimes

A literal key-set audit partitions the post-braid set without overlap:

```text
global weighted AM-GM certificates        1,623
former fan keys, now global AM-GM            248
translated-orthant proof trees               387
uniform-margin finite triangles                2
existing structural theorems                   6
                                            -----
total                                       2,266
```

### Weighted AM-GM

For a residual exponential sum, each negative term receives rational weights
on positive spectral terms. The replay checks

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l),
sum_n alpha_(n,j) D_n <= C_j.
```

Capacities use exact rational intervals. Products are checked after taking
logs with 512-bit directed rounding. The two ledgers replay as

```text
O11_AMGM_EXACT PASS records=1623 threads=5
O11_AMGM_EXACT PASS records=248 threads=5
```

### Translated orthants

The remaining rectangular certificates recursively use the exact identity

```text
N^d = {coordinate=0 face} disjoint union {coordinate>=1 shifted tail}.
```

Every terminal node is either an exact AM-GM certificate or an integral fusion
leaf. The complete forest reports

```text
O11_TREE_EXACT PASS
roots=387 amgm=1394 splits=1207 leaves=200 threads=5
```

### Uniform-margin triangles

Two two-variable chambers admit a common log-base margin `delta=0.14`.
Outside finite triangles, exponential slope dominates the total negative
coefficient mass. The remaining lattice points are evaluated by exact fusion
dynamic programming:

```text
O11_MARGIN2_EXACT PASS
cases=2 delta=0.14 finite_checks=7068 finite_min=48 threads=5
```

### Structural regimes

Two keys belong to the already-published first `O_11` chamber theorem. Four
keys are instances of the all-rank separated-pole theorem. The key partition
replay reports

```text
O11_KEY_PARTITION_EXACT PASS
post_braid=2266 amgm=1623 former_fan=248 tree=387 margin=2 structural=6
```

## 4. Theorem

**Theorem 4.1.** The rank-six orbit ring `O_11` satisfies scalar `GKS2*` for
signed words of arbitrary length.

**Proof.** Support overlap reduces every word to the support-disjoint parity
chambers enumerated above. The exact universal census and post-braid key audit
cover every one of the 27,962 residual regimes. Each regime is certified by a
nonnegative spectral sum, an exact transport, weighted AM-GM, a translated
orthant tree, a uniform-margin finite check, or an existing structural theorem.
Thus every even-minus scalar corner is nonnegative. Reversing support overlap
proves the result for all words. QED.

**Corollary 4.2.** The complete fusion ring `SU(2)_11` satisfies scalar
`GKS2*` for every signed word of arbitrary length. Every partial-character
coefficient is nonnegative.

**Proof.** Apply the exact odd-level simple-current lift, followed by the
plus-factor identity for the partial-character column. QED.

## 5. Ordinary `SU(2)` stable consequence

For an ordinary signed word with total label sum `T` and target `V_c`, put

```text
K(c)=max(max_i a_i,c,ceil((T+c)/2)).
```

The finite theorems now prove every ordinary partial coefficient satisfying

```text
K(c)<=11.
```

For the Ginibre corner this becomes

```text
max(max_i a_i,ceil(T/2))<=11.
```

In particular, every ordinary word with labels at most eleven and total label
sum at most twenty-two is covered, together with the additional words allowed
by the sharper combined bound.

## 6. Reproduction

Run

```text
character_ring_iter/replay_su2_o11_full_gks2.sh
```

from a checkout of this repository. The wrapper reconstructs the deterministic
certificate archive, verifies its SHA-256, replays every new verifier, and also
runs the already-published first-chamber verifier.

The next unresolved odd-level finite target is the rank-seven orbit ring
`O_13`, equivalent by simple-current lifting to full `SU(2)_13`.
