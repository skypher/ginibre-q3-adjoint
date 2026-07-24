# SU(2) finite-fusion reduction and exact stable passage

Date: 2026-07-24

## Status

This note proves three reductions for the character-cone problem:

1. labels occurring with both signs can be eliminated inductively;
2. every irreducible difference factors through the fundamental-character
   difference with a coefficientwise nonnegative kernel; and
3. every fixed signed word in the ordinary SU(2) representation ring is
   reproduced exactly in SU(2)_k once the level is at least the sum of its
   highest-weight labels.

The third point makes the passage from a uniform finite-level theorem to the
ordinary SU(2) theorem completely formal.  This note does **not** prove the
remaining uniform finite-level `GKS2*` inequality.

## 1. Signed coefficients

Write `V_a` for the irreducible SU(2) representation of highest weight `a`
and `chi_a` for its character.  For labels `a_1,...,a_r` and signs
`epsilon_i in {+1,-1}`, put

```text
C_infinity(a,epsilon;c)
 = [V_c tensor V_0]
   product_i (V_(a_i) tensor V_0
              + epsilon_i V_0 tensor V_(a_i))
```

in `R(SU(2)) tensor R(SU(2))`.  Its corner `c=0` is the signed Haar integral
appearing in Ginibre's `Q_3`; nonnegativity for every `c` is the stronger
partial-character form used in the finite-level program.

At level `k`, use the same definition in the fusion ring `R_k`, whose product
is

```text
V_a *_k V_b
 = sum V_c,
   |a-b| <= c <= min(a+b,2k-a-b),
   c == a+b (mod 2).
```

Denote the corresponding coefficient by `C_k(a,epsilon;c)`.

## 2. Elimination of overlapping sign supports

Suppose the same label `a` occurs once with each sign.  In either the
ordinary representation ring or a finite fusion ring,

```text
(V_a tensor 1-1 tensor V_a)(V_a tensor 1+1 tensor V_a)
 = V_a^2 tensor 1-1 tensor V_a^2.
```

The square `V_a^2` is a nonnegative integral combination of simple objects.
Distributing this combination expresses the original signed coefficient as
a nonnegative sum of signed coefficients with one fewer factor.  Induction
therefore reduces `Q_3` and `GKS2*` to words in which no irreducible label
occurs with both signs.

This is the support-disjoint chamber used throughout the computational and
structural reductions in this repository.

## 3. Fundamental-difference factorization

Let `x=chi_1(g)` and `y=chi_1(h)`.  Since

```text
chi_a = U_a(chi_1/2),
```

where `U_a` is the Chebyshev polynomial of the second kind, the elementary
divided-difference identity gives

```text
chi_a(g)-chi_a(h)
 = (chi_1(g)-chi_1(h))
   sum_(j=0)^(a-1) chi_j(g) chi_(a-1-j)(h).          (3.1)
```

Indeed, both sides satisfy the same three-term recurrence in `a`, and the
identity is immediate for `a=1,2`.  Equivalently it is the standard formula

```text
(U_a(X)-U_a(Y))/(2X-2Y)
 = sum_(j=0)^(a-1) U_j(X)U_(a-1-j)(Y).
```

Every coefficient in the kernel in (3.1) is nonnegative.  Hence a word with
`2m` minus factors has the exact form

```text
(chi_1(g)-chi_1(h))^(2m)
 product_(i in minus) K_(a_i)(g,h)
 product_(j in plus) (chi_(a_j)(g)+chi_(a_j)(h)),
```

where

```text
K_a(g,h)=sum_(j=0)^(a-1) chi_j(g)chi_(a-1-j)(h).
```

This isolates the only sign-changing factor.  It is a reduction, not by
itself a positivity proof, because the remaining kernels couple the two
variables asymmetrically.

## 4. Exact stabilization theorem

**Theorem.**  Let `a_1,...,a_r` be nonnegative labels and set

```text
T=a_1+...+a_r.
```

For every sign pattern, every target `c`, and every level `k>=T`,

```text
C_k(a,epsilon;c)=C_infinity(a,epsilon;c).            (4.1)
```

**Proof.**  Expand the signed product into ordinary products in the two
tensor factors.  It is enough to compare one product of a submultiset of the
labels in `R_k` with its ordinary Clebsch--Gordan product.

Choose any binary parenthesization.  Attach to each node the sum of the leaf
labels below that node.  Inductively, every simple label occurring at that
node is at most that attached sum.  If two child nodes have attached sums
`A,B` and carry simple labels `u<=A`, `v<=B`, then

```text
u+v <= A+B <= T <= k.
```

Consequently

```text
2k-u-v >= u+v,
```

so the affine upper cutoff in the fusion rule is inactive:

```text
min(u+v,2k-u-v)=u+v.
```

Thus that multiplication step is exactly the ordinary Clebsch--Gordan step.
Induction up the tree proves equality of the complete products.  Applying
this independently in the two tensor factors and then summing the signed
expansion proves (4.1).  QED.

## 5. Cofinal finite-to-stable implication

**Corollary.**  Assume `GKS2*` holds in `SU(2)_k` for every finite level `k`:
for every signed word with an even number of minus signs, its corner
coefficient is nonnegative.  Then Ginibre's `Q_3` holds for the cone generated
by all irreducible real SU(2) characters.

More strongly, if every finite fusion ring has nonnegative partial-character
coefficients `C_k(a,epsilon;c)`, then the same partial-character statement
holds in the ordinary representation ring.

**Proof.**  For a fixed word choose any `k>=T`.  The desired ordinary
coefficient equals the finite-level coefficient by the theorem and is
therefore nonnegative.  Ginibre's generator-to-cone argument then promotes
the inequalities from the irreducible generating family to its closed
multiplicative convex cone.  QED.

The limit is therefore not analytic: each fixed coefficient is eventually
constant in `k`.

## 6. Remaining finite target

After support-disjoint reduction, the unresolved statement can be written
purely in the based ring `R_k`:

```text
[V_0 tensor V_0]
 product_i (V_(a_i) tensor V_0
            + epsilon_i V_0 tensor V_(a_i)) >= 0,     (GKS2*)
```

for every level `k`, every word, and every even minus parity.  The stronger
surviving target is

```text
[V_c tensor V_0](same product) >= 0                  (PC)
```

for every `0<=c<=k`.

The repository already proves substantial sectors, including the universal
low-factor range, the complete ordinary six-factor sector, all-length
results for the first odd simple-current orbit rings, and the uniform
top-orbit ray.  It also records exact counterexamples to two stronger but
unnecessary approaches: gamma-positivity and log-concavity of the associated
palindromic coefficient polynomial.  Thus the next proof should target
`(GKS2*)` or `(PC)` directly, most plausibly through the cyclic shell/outer-
Turan transport developed in `CENTRAL_CHARACTER_Q3_SEARCH.md`.

## 7. Bounded regression

`character_ring_iter/verify_su2_finite_stable_transfer.py` independently
checks (4.1) by exact integer dynamic programming in a finite box.  The
regression is only a check of the implementation and indexing; the theorem
above is the unbounded proof.
