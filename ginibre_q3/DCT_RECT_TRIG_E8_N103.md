# DCT-RectTrig(E8, n=103) certificate

This note records the wide-cap rectangular certificate for the single E8
bridge step

```text
D_E8(103) = Q_3^{E8}(105) - 4 Q_3^{E8}(103).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities and cubic scalar trigonometric bounds are those of
Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker `direct_chain_rect_e8_n103_certificate.py` imports the
same exact Sturm and root-identity checks and returns

```text
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
root_identity_constants: OK
```

## Lemma 2 (quartic Weyl-sine and exponential bounds)

The quartic Weyl-sine bound

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100,     |z| <= 4,
```

and the fractional-part exponential lower bound using `e < 1457/536` are
those of Lemma 2 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker returns

```text
sine_quartic_sturm: OK
e_upper_1457_536: OK
```

## Lemma 3 (shifted Chebyshev degree-14 negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character, set
`S = X+Y`, and set

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

Let

```text
C_14(S) = T_14((S-249)/247) / T_14(-265/247),
```

where `T_14` is the Chebyshev polynomial of the first kind.  Then
`C_14(-16)=1`.  On the negative-sign region for odd `n=103`,

```text
S in [-16,-2] union [0,2].
```

For every such `S`,

```text
|S^103(S^2-4)|
 <= (63/65) * 16^95 * S^8 * (S^2+4) * C_14(S)^2.
```

Proof.  The right side is non-negative for all real `S`.  The affine
argument `(S-249)/247` maps `[2,496]` to `[-1,1]`, so no zero of `T_14`
lies in the positive-sign interval `[0,2]`.  On `S in [-16,-2]`, put
`x=-S`.  After clearing positive denominators and dividing by
`x^8(16-x)`, the checker verifies by exact Sturm root count that the
residual quotient is positive on `[2,16]`.  On `S in [0,2]`, after
clearing positive denominators and dividing by `S^8`, the checker verifies
by exact Sturm root count that the residual quotient is positive on
`[0,2]`.  It returns

```text
shifted_chebyshev_normalization: OK
shifted_chebyshev_negative_interval_sturm: OK
shifted_chebyshev_positive_interval_sturm: OK
```

The exact moment prefix `m_0,...,m_40` from
`character_ring_iter/logs/e8_70.log`, with OEIS overlap `m_0,...,m_30`
checked against `references/oeis_A179663.txt`, gives

```text
K_38 =
  int int (X-Y)^2 S^8(S^2+4) C_14(S)^2 dX dY
  = 211367009391307597006963601980640081834787585716146581122825466841199954628
    / 101897318230030856294954778256196368953775863550819861054942528948894561.
```

Therefore the negative part of `D_E8(103)` is bounded by

```text
Neg(103) = (63/65) * 16^95 * K_38.
```

The checker returns

```text
moment_log_0_40_oeis_overlap_0_30: OK
shifted_chebyshev_moment_integral: OK
```

## Lemma 4 (wide-cap rectangular union)

Let `U_103` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 199/4},
t0 in {30, 601/20, ..., 959/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`56387` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 103 / (4*248) = 1545/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n103_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_103`,

```text
chi(y) - chi(x) >= (496/103) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/103
       + L_3 (496^2/1800)s0^3/103^2
  - U_1 t1 - U_2 (496/50)t1^2/103
       - U_3 (496^2/1800)t1^3/103^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/103
  + C_2(496/50)(s0^2+t0^2)/103^2
  + C_3(496^2/1800)(s1^3+t1^3)/103^3.
```

Moreover every retained cell has `g_R > 0`, `h_R > 0`, and
`496^2 h_R^2 - 4 > 0`.

Proof.  Lemma 1 gives the displayed formulas.  The endpoint choices use
lower endpoints in positive terms and upper endpoints in negative terms.
The exact rational checker verifies the positivity comparisons over all
retained cells.

## Proposition 6 (rectangular radial mass lower)

For a retained cell `R=[s0,s1] x [t0,t1]`, define

```text
F_R =
  g_R^2 * ExpLower(s1+t1)
  * (s1^124-s0^124)/124
  * (t1^124-t0^124)/124
  * 1/(123!)^2
  / 124.
```

This is a lower bound for the normalized radial mass factor on `R`.

Proof.  On the cell the gap factor is at least `g_R^2`; Lemma 2 supplies
the radial Gaussian lower bound.  Since `a=124`, `Gamma(124)=123!`, giving
the displayed factor.

## Lemma 7 (exact-rational arithmetic check)

Let

```text
A0_lower = A0_prefactor(E8) * (113/355)^8,
Sine_R =
  ExpLower(
    2(1/24) * 496(s1+t1)/103
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/103^2
  ),
Neg(103) = (63/65) * 16^95 * K_38.
```

The rational lower ratio

```text
R103 =
  sum_R A0_lower * 496^103 * 103^(-250)
        * F_R * h_R^103
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(103)
```

satisfies `R103 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n103_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=103 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=50 t_start=30 t_end=48
retained_cells=56387 cap=1545/31
shifted_chebyshev_degree=14 cheb_at_0=9596503446458386511349630846764273/319213593429275590263910293224524081
cheb_integral=211367009391307597006963601980640081834787585716146581122825466841199954628/101897318230030856294954778256196368953775863550819861054942528948894561
log10_cheb_negative_vs_p12=-3.92
negative_bound=(63/65)*16^95*cheb_integral
log10_union_fraction_lower=-39.97
gap103_min=2702191337/19096200000000
h103_min=3709027066736679977/10585792812500000000 direct_factor103_min=13755059802573067332697806102512720529/455494803058242797851562500000000
sine_exponent_max=47992142468/1160359375
log10_ratio103_lower=0.68
log10_best_cell_ratio103=-3.31
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
sine_quartic_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
moment_log_0_40_oeis_overlap_0_30: OK
shifted_chebyshev_normalization: OK
shifted_chebyshev_negative_interval_sturm: OK
shifted_chebyshev_positive_interval_sturm: OK
shifted_chebyshev_moment_integral: OK
cell_count_56387: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h103: OK
positive_direct103: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio103_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=103))

The rectangular-union positive contribution to `D_E8(103)` dominates the
negative-region bound `Neg(103)`.  Thus `D_E8(103) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=103`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that
the lower bound is larger than the negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=103)

Assume `E8-FiniteBridge(103)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 103`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(105) >= 4 Q_3^{E8}(103)`, hence the
assumed finite bridge through `n=103` gives `Q_3^{E8}(105) >= 0`.  This
supplies the hypothesis `E8-FiniteBridge(105)`, and Corollary 9 of
`DCT_RECT_TRIG_E8_N105.md` applies.
