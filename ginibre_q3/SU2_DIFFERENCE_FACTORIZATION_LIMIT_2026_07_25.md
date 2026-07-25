# The difference-factorization reduction has exactly one ordinary instance

Date: 2026-07-25

## Motivation

Two results dated 2026-07-25 use the same algebraic move.

`SU2_V2_MINUS_FUNDAMENTAL_REDUCTION_2026_07_25.md` uses

```text
V_1^2=V_0+V_2                =>  D_2=D_1 S_1,
```

and `SU2_ODD_ORBIT_B1_TOP_REDUCTION_2026_07_25.md` uses, for the top orbit
label `T=B_(m-1)` of `O_(2m-1)`,

```text
T^2=B_0+B_1                  =>  D_1=D_T S_T.
```

In both cases a minus factor becomes one minus and one plus factor of another
label. Minus parity is preserved, and the reduction is exact: no spectral
estimate, no chamber census, and uniform in both level and word length. Those
are exactly the properties the level-by-level program does not have, so it is
worth knowing how far the move generalises.

It does not generalise. In the ordinary representation ring it has exactly one
instance, and the orbit-ring instance is a truncation phenomenon.

## 1. Statement

Write `D_a=chi_a(x)-chi_a(y)` and `S_a=chi_a(x)+chi_a(y)`.

**Theorem.** In the ordinary `SU(2)` representation ring,

```text
D_a=D_u S_u
```

holds if and only if `(a,u)=(2,1)`.

**Proof.** `D_u S_u=chi_u(x)^2-chi_u(y)^2`, so the identity holds for all
`x,y` exactly when

```text
chi_u(x)^2-chi_a(x)=chi_u(y)^2-chi_a(y)
```

for all `x,y`, that is, when `chi_u^2-chi_a=c` for a constant `c`.

Clebsch--Gordan gives

```text
chi_u^2=sum_(i=0)^(u) chi_(2i),                        (1.1)
```

a sum of exactly `u+1` distinct irreducible characters. The only constant
irreducible character is `chi_0`. So `chi_u^2-chi_a` is constant only if the
sum (1.1) has at most two terms, one of which is `chi_a` and the other
`chi_0`. Hence `u+1<=2`. The case `u=0` is trivial, so `u=1`, and then

```text
chi_1^2=chi_2+chi_0
```

gives `a=2` and `c=1`. QED

## 2. Why the orbit ring admits a second instance

The proof uses (1.1), which is the *untruncated* Clebsch--Gordan rule. In a
fusion ring the square of the top label truncates to two terms, which is
precisely the hypothesis the theorem otherwise forbids. So the orbit-ring
reduction is available once per level, at the top label only, and is not an
additional ordinary-ring identity.

## 3. The weaker multiplicative question

Since `D_a=D_1 K_a` with the coefficientwise nonnegative kernel

```text
K_a=sum_(j=0)^(a-1) chi_j(x) chi_(a-1-j)(y),
```

one may ask the weaker question of whether `K_a` factors into a product of
`S_b`, which would still give a parity-preserving reduction. A symbolic
factorisation over `Q[x,y]` for `2<=a<=16` finds

```text
a=2:       K_2=S_1;
a even:    K_a=S_1 * (irreducible), the cofactor never an S_b;
a odd:     K_a irreducible.
```

For example `K_4=S_1(x^2+y^2-3)` while `S_2=x^2+y^2-2`. So `a=2` is the only
case in that range here too. This search is supporting evidence only; the
theorem of Section 1 is the unbounded statement.

## 4. Consequence

The `{1,2}` label sector proved in
`SU2_V1_V2_SUPPORT_DISJOINT_GKS_2026_07_25.md` is maximal for this technique.
There is no bootstrap from label `2` to label `3` by the same factorisation,
so the label axis does not extend by this route and effort is better spent on
the factor axis or on the Walsh-coefficient target.

This is a negative result. Its value is that it closes a direction cheaply.

## 5. Exact verification

`character_ring_iter/verify_su2_difference_factorization_limit.py` checks,
in exact integer bivariate polynomial arithmetic with no external dependency:

1. `chi_u^2=sum_(i=0..u) chi_2i` for `u<=24`;
2. `D_2=D_1 S_1` exactly;
3. `D_a != D_u S_u` for every other pair with `a,u<=48`.

It reports

```text
SU2_DIFFERENCE_FACTORIZATION_LIMIT PASS
clebsch_gordan=25 pairs=2303 instances=1 unique_pair=(a=2,u=1).
```

The bounded run is a cross-check; Section 1 is the proof.
