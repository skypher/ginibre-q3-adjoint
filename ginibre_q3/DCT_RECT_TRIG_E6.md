# DCT-RectTrig(E6) rank-corrected certificate

This note records the active exact certificate for

```text
D_E6(n) = Q_3^E6(n+2) - 4 Q_3^E6(n),   odd n >= 75.
```

It supersedes the former single-rectangle note, which used the false bound
`chi_ad >= -2`.  Garibaldi--Guralnick--Rains, Theorem 2.1, gives the valid
uniform bound `chi_ad >= -rank(E6) = -6`; their Theorem 2.3 and Table 1 give
the sharper actual minimum `-3`.  The active proof deliberately uses the
uniform rank bound, so its negative-region estimate is

```text
U_E6(n) = 2*((12)^2+4)*12^n = 296*12^n.
```

## Lemma 1 (root and scalar identities)

For E6,

```text
d=78, r=6, |R_+|=36, kappa=12, alpha=80,
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/16.
```

The quartic identity follows because the basic Weyl-invariant degrees are
`2,5,6,8,9,12`, so every invariant quartic is a scalar multiple of `Q^2`;
evaluation at a root gives the denominator `16`.

For the adjoint global form, the normalized two-torus Weyl coefficient has
the additional lattice factor

```text
nu_E6 = |P^vee/Q^vee|^2 / covol(Q^vee)^2 = 3^2/3 = 3.
```

Here both the Cartan determinant and the squared coroot-lattice covolume are
`3`.  The checker verifies the determinant before consuming this factor.

## Lemma 2 (wide-cap scalar bounds)

On `|z| <= 14/3`, the exact Bernstein checks in
`direct_chain_rect_e6_e7_rank_certificate.cpp` prove

```text
2(1-cos z) >= 0.99119225 z^2 - 0.07802858 z^4
                  + 0.00169078 z^6,
2 cos z >= 2-z^2+(83/2000)z^4,
sin(z/2)/(z/2) >= exp(-z^2/24-z^4/2880-z^6/103800).
```

The checker compares the rational polynomials with alternating Taylor
minorants on four exact subintervals.  No floating-point sign is used.

## Proposition 3 (OpenMP/GMP rectangle union)

At `n0=75`, partition

```text
10 <= s < 30,   5 <= t < 20
```

into cells of side `1/10`.  Retain exactly the `23088` cells on which the
certified character gap and character base are positive, both are
nondecreasing for all later `n`, and the exact odd-step ratio is greater
than one.  Every cell lies in the root-angle cap because

```text
s,t < (14/3)^2*12*75/(4*78) = 2450/39.
```

The checker integrates rational lower endpoint densities with the exact
integer Gamma factor, the Macdonald--Mehta normalization, and rational
lower bounds for both exponentials.  Its accepted replay gives

```text
retained_cells=23088
log10_positive_lower=92.314
log10_negative_upper=83.410
log10_base_ratio_lower=8.904
log10_odd_step_ratio_lower=0.573
case_result: PASS
```

Consequently `D_E6(n) >= 0` for every odd `n >= 75`.  The exact prefix has
`Q_3^E6(75)>0`, so Chain propagation proves the complete odd tail.
