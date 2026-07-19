# DCT-RectTrig(E7) rank-corrected certificate

This note records the active exact certificate for

```text
D_E7(n) = Q_3^E7(n+2) - 4 Q_3^E7(n),   odd n >= 65.
```

It supersedes the former single-rectangle note, which used the false bound
`chi_ad >= -3`.  Garibaldi--Guralnick--Rains, Theorem 2.1, gives the valid
bound `chi_ad >= -rank(E7) = -7`.

## Lemma 1 (root and scalar identities)

For E7,

```text
d=133, r=7, |R_+|=63, kappa=18, alpha=135,
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/27.
```

The invariant-degree argument and evaluation at a root give the quartic
denominator `27`.  The wide-cap scalar inequalities are the exact rational
ones recorded in Lemma 2 of `DCT_RECT_TRIG_E6.md`.

For the adjoint global form, normalized Haar measure contributes the exact
two-torus lattice factor

```text
nu_E7 = |P^vee/Q^vee|^2 / covol(Q^vee)^2 = 2^2/2 = 2.
```

The exact checker verifies the Cartan determinant `2` before using it.

## Lemma 2 (degree-26 negative-region majorant)

Let `psi_E7` be the Chain pushforward measure, `S` its coordinate, and

```text
P(S) = T_26((2S-273)/259) / T_26((-28-273)/259).
```

At `n0=65`, exact Bernstein coefficients prove on the two negative-sign
intervals `S in [-14,-2]` and `S in [0,2]` that

```text
|S^65(S^2-4)|
 <= (24/25)*14^(65-14)*S^14*(S^2+4)*P(S)^2.
```

The right side integrates exactly from the Chain moments.  Its degree is at
most `68`, so it uses only the exact E7 adjoint prefix `m_0,...,m_70`.
Multiplication by `14^2` propagates the majorant from `n` to `n+2` because
`|S|<=14` on both sign intervals.

## Proposition 3 (OpenMP/GMP rectangle union)

At `n0=65`, partition

```text
15 <= s < 35,   8 <= t < 28
```

into cells of side `1/20`.  Retain exactly the `66925` cells on which the
gap and character base are nondecreasing and the exact positive-side
odd-step ratio exceeds `101/100`.  Each cell lies in the cap
`s,t < (14/3)^2*18*65/(4*133) = 910/19`.

The C++ checker uses only GMP integers and rationals.  For the half-integer
radial power it bounds each square root downward by an exactly verified
six-decimal rational.  The accepted replay gives

```text
retained_cells=66925
log10_positive_lower=61.666
log10_negative_upper=61.240
log10_base_ratio_lower=0.426
odd_step_ratio > 101/100
case_result: PASS
```

Therefore `D_E7(n) >= 0` for every odd `n >= 65`.  The exact prefix has
`Q_3^E7(65)>0`, so Chain propagation proves the complete odd tail.
