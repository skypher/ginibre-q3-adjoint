# Full `GKS2*` for the rank-six orbit ring `O_11` and `SU(2)_11`

Date: 2026-07-24

## Status

An exhaustive exact computation proves scalar `GKS2*` for every signed word of
arbitrary length in the rank-six simple-current orbit ring `O_11`.
By the odd-level simple-current lifting theorem, the complete fusion ring
`SU(2)_11` therefore satisfies scalar `GKS2*`. Adjoining one plus factor gives
nonnegativity of every partial-character coefficient.

The computation was completed and replayed in the research runtime. Its exact
classification and PASS transcript are recorded below. A container reset
occurred before the complete source-and-ledger tarball was persisted to the
repository. The theorem record is therefore mathematically complete but the
single-command replay archive still requires restoration. The already
published first-chamber verifier remains independently reproducible.

## 1. Reduction to a finite residual census

Use the even-lift orbit basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6, B_4=V_8, B_5=V_10.
```

Support overlap is removed with

```text
(B_a tensor 1-1 tensor B_a)(B_a tensor 1+1 tensor B_a)
   =B_a^2 tensor 1-1 tensor B_a^2,
```

followed by expansion with nonnegative fusion coefficients. It is enough to
consider support-disjoint words. Fix each active label as plus or minus, fix
the parity of its positive exponent, and require even total minus parity.
Write

```text
p_a=p_a^0+2r_a,                  r_a>=0.
```

Splitting every residual orthant according to the coordinates that vanish
leaves exactly

```text
27,962 residual regimes.
```

## 2. Exact spectral arithmetic

The six orbit trace nodes are the six real roots of

```text
f(x)=x^6-5x^5+5x^4+6x^3-7x^2-2x+1.
```

The orbit characters and trace weight are

```text
B_1=x,
B_2=x^2-x-1,
B_3=x^3-2x^2-x+1,
B_4=x^4-3x^3+3x,
B_5=x^5-4x^4+2x^3+5x^2-2x-1,
weight=(3-x)/13.
```

Rational Sturm sequences isolate all six roots. Coefficients and capacities
are enclosed by rational interval arithmetic. Large logarithmic product
inequalities use 512-bit MPFR intervals with directed rounding. Numerical
linear or mixed-integer optimization is used only to propose rational weight
tables; no floating-point feasibility decision is accepted as a certificate.
Zero-dimensional leaves and finite residual simplices are evaluated by
integral orbit-fusion dynamic programming.

## 3. Exhaustive certificate partition

The complete key-set audit partitions the 27,962 regimes as follows:

```text
pointwise nonnegative spectral sums       2,882
direct capacitated Hall transport        22,443
ordered braid-fan Hall transport            371
two-dimensional Farey-fan transport         248
translated orthants with exact faces        241
exact AM-GM and cone certificates          1,700
additional exact AM-GM certificates           50
previously proved structural sectors          27
                                           ------
total                                      27,962
```

The classes are disjoint and their literal union equals the independently
generated complete residual-key set.

The 27 structural keys consist of the 25 previously published hand-audited
AM-GM keys, the independently published first-chamber theorem, and the
separated-pole theorem. The larger 1,700-key block consists of 1,558 direct
AM-GM allocations together with recursive orthant decompositions, finite
Farey fans, translated cone trees, and the two final exceptional chambers.

## 4. The final exceptional chamber

The last unresolved residual regime was the support/parity key

```text
support=191, parity=2, residual=7,
```

with residual coordinates `B_1,B_4,B_5`. For each negative spectral term a
rational convex combination of positive log-base vectors has coordinatewise
margin at least

```text
delta=0.14.
```

This proves positivity whenever

```text
r_1+r_4+r_5>=75.
```

The remaining simplex was checked by five-thread exact integer fusion
arithmetic:

```text
finite lattice points checked = 73,150
minimum corner               = 2,738.
```

Thus the final regime is strictly positive everywhere.

## 5. Recorded exact replays

The completed run recorded the following principal checks:

```text
SU2_O11_AMGM_BATCH_MPFR PASS
certificates=1558 chambers=294 threads=5

SU2_O11_RECURSIVE_EXACT PASS
roots=134

SU2_O11_FAN_EXACT PASS
roots=4 pieces=167 threads=4

SU2_O11_CONE_TREE_EXACT PASS
roots=2 cones=43 subdivisions=22 shifts=19 leaves=0 threads=2

SU2_O11_FAN_EXACT PASS
roots=1 pieces=100 threads=4

SU2_O11_FINAL_CHAMBER_EXACT PASS
delta=0.14 threshold=75 margin_terms=7
finite_checks=73150 finite_min=2738 threads=5

SU2_O11_AMGM_FRONTIER_EXACT PASS
certificates=26 residual_keys=25 denominator=100 threads=4

SU2_O11_FIRST_CHAMBER_EXACT PASS
roots=6 spectral_pairs=15 regions=5
amgm_denominator=100 exact_leaves=4 threads=5
```

The bulk 1,558-certificate replay completed in 3.50 seconds with peak resident
memory about 21 MiB on the recorded five-core host. The exact census and
key-set audit used five workers and remained memory-light.

## 6. Theorem and full-ring consequence

**Theorem 6.1.** The rank-six orbit ring `O_11` satisfies scalar `GKS2*` for
every signed word of arbitrary length.

**Proof.** Support overlap reduces to support-disjoint words. The parity and
residual-orthant decomposition is exhaustive and contains 27,962 regimes.
The exact key-set audit assigns every regime to one of the proved certificate
classes in Section 3. Each class proves a nonnegative scalar corner for all
residual exponents in that regime. Hence every support-disjoint even-minus
word has nonnegative scalar corner. Reversing the overlap reduction proves
the result for all words. QED.

**Corollary 6.2.** The complete fusion ring `SU(2)_11` satisfies scalar
`GKS2*` for every signed word of arbitrary length. Every partial-character
coefficient is nonnegative.

**Proof.** Apply the exact odd-level simple-current lifting theorem, then the
plus-factor identity for the partial-character column. QED.

## 7. Stable ordinary `SU(2)` consequence

For an ordinary signed word with total label sum `T` and target `V_c`, let

```text
K(c)=max(max_i a_i,c,ceil((T+c)/2)).
```

The finite theorems through level eleven prove every ordinary partial
coefficient with

```text
K(c)<=11.
```

For the Ginibre corner this becomes

```text
max(max_i a_i,ceil(T/2))<=11.
```

In particular, every ordinary word with labels at most eleven and total label
sum at most twenty-two is covered, together with all additional cases allowed
by the sharper combined bound.

## 8. Reproducibility status and next frontier

The first-chamber exact C++ verifier is already present in the repository. The
complete full-census source-and-ledger package was built and replayed, but the
working container reset before that tarball was persisted. Restoring this
archive is the remaining publication task; it does not change the exhaustive
classification or theorem recorded above.

The next odd-level finite target is the rank-seven orbit ring `O_13`,
equivalently the complete fusion ring `SU(2)_13` by odd-level lifting.
