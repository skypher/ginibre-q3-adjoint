# DCT-RectTrig(E8, n=99 and n=97) certificate

This note records the wide-cap rectangular certificate for the two E8
bridge steps

```text
D_E8(99) = Q_3^{E8}(101) - 4 Q_3^{E8}(99),
D_E8(97) = Q_3^{E8}(99) - 4 Q_3^{E8}(97).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities and cubic scalar trigonometric bounds are those of
Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker `direct_chain_rect_e8_n99_n97_certificate.py` imports
the same exact Sturm and root-identity checks and returns

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

## Lemma 3 (shifted Chebyshev degree-26 negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character, set
`S = X+Y`, and set

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

Let

```text
C_26(S) = T_26((S-249)/247) / T_26(-265/247),
```

where `T_26` is the Chebyshev polynomial of the first kind.  Then
`C_26(-16)=1`.  For each `n in {99,97}`, on the negative-sign region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^n(S^2-4)|
 <= (63/65) * 16^(n-8) * S^8 * (S^2+4) * C_26(S)^2.
```

Proof.  The right side is non-negative for all real `S`.  On
`S in [-16,-2]`, put `x=-S`.  After clearing positive denominators and
dividing by `x^8(16-x)`, the checker verifies positivity on `[2,16]` by
checking that every Bernstein coefficient of the residual quotient is
positive.  On `S in [0,2]`, after clearing positive denominators and
dividing by `S^8`, the checker verifies positivity on `[0,2]` by the same
Bernstein-coefficient test.  It returns

```text
shifted_chebyshev_normalization: OK
shifted_chebyshev_negative_bernstein_99: OK
shifted_chebyshev_positive_bernstein_99: OK
shifted_chebyshev_negative_bernstein_97: OK
shifted_chebyshev_positive_bernstein_97: OK
```

The exact moment prefix `m_0,...,m_64` from
`character_ring_iter/logs/e8_70.log`, with OEIS overlap `m_0,...,m_30`
checked against `references/oeis_A179663.txt`, gives

```text
K_62 =
  int int (X-Y)^2 S^8(S^2+4) C_26(S)^2 dX dY
  = 652708718049108501279891206610792628257515262318580914182134886342747101135881782357617401218658855060065796860047212614264650354628
    / 2445200692443756455056633566892545977519407487004403208271763856863104974462205280301950124320130386149081237313785855952197768855841.
```

Therefore the negative part of `D_E8(n)` is bounded by

```text
Neg(n) = (63/65) * 16^(n-8) * K_62,     n in {99,97}.
```

The checker returns

```text
moment_log_0_64_oeis_overlap_0_30: OK
shifted_chebyshev_moment_integral: OK
```

## Lemma 4 (wide-cap rectangular unions)

For `n=99`, let `U_99` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 957/20},
t0 in {30, 601/20, ..., 957/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`42662` cells.

For `n=97`, define `U_97` the same way, with endpoint sets

```text
s0 in {40, 801/20, ..., 937/20},
t0 in {30, 601/20, ..., 937/20}.
```

The exact enumeration has `35869` cells.  Each retained cell lies inside the
root cap `|alpha.v| <= 4`, whose scalar cap is

```text
15n/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n99_n97_certificate.py`.

## Lemma 5 (local character lower bounds on the unions)

For `n in {99,97}` and for a retained cell `R=[s0,s1] x [t0,t1]` in
`U_n`,

```text
chi(y) - chi(x) >= (496/n) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/n
       + L_3 (496^2/1800)s0^3/n^2
  - U_1 t1 - U_2 (496/50)t1^2/n
       - U_3 (496^2/1800)t1^3/n^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/n
  + C_2(496/50)(s0^2+t0^2)/n^2
  + C_3(496^2/1800)(s1^3+t1^3)/n^3.
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
Sine_R(n) =
  ExpLower(
    2(1/24) * 496(s1+t1)/n
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/n^2
  ),
Neg(n) = (63/65) * 16^(n-8) * K_62.
```

For `n in {99,97}`, the rational lower ratio

```text
R_n =
  sum_{R in U_n} A0_lower * 496^n * n^(-250)
        * F_R * h_R^n
        * (496^2 h_R^2 - 4)
        * Sine_R(n)
  / Neg(n)
```

satisfies `R_n > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n99_n97_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n_values=99,97 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 t_start=30 shifted_chebyshev_degree=26 moment_max=64
cheb_at_0=2216397178979222624328797476162428806425206380934057202670732753/1563713750161376708182872704466724144181174906564989485936798408721
cheb_integral=652708718049108501279891206610792628257515262318580914182134886342747101135881782357617401218658855060065796860047212614264650354628/2445200692443756455056633566892545977519407487004403208271763856863104974462205280301950124320130386149081237313785855952197768855841
case n=99 s_end=479/10 t_end=479/10
retained_cells=42662 cap=1485/31
negative_bound=(63/65)*16^91*cheb_integral
log10_cheb_negative_vs_p12=-7.81
log10_union_fraction_lower=-42.52
gap99_min=1856998249/13750000000000
h99_min=3540441037296109973/10254296250000000000 direct_factor99_min=12533013084050604192273505118110060729/427413629937744140625000000000000
sine_exponent99_max=134717388541/3215953125
log10_ratio99_lower=1.70
log10_best_cell_ratio99=-2.15
case n=97 s_end=469/10 t_end=469/10
retained_cells=35869 cap=1455/31
negative_bound=(63/65)*16^89*cheb_integral
log10_cheb_negative_vs_p12=-7.81
log10_union_fraction_lower=-43.97
gap97_min=11851190925583/105851250000000000
h97_min=47220675071042327/136724531250000000 direct_factor97_min=2229488212999415671529399725574929/75985291385650634765625000000
sine_exponent97_max=25857736616/617465625
log10_ratio97_lower=0.27
log10_best_cell_ratio97=-3.50
moment_log_0_64_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein_99: OK
shifted_chebyshev_positive_bernstein_99: OK
ratio99_gt_1: OK
shifted_chebyshev_negative_bernstein_97: OK
shifted_chebyshev_positive_bernstein_97: OK
ratio97_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=99 and n=97))

The rectangular-union positive contribution to each of `D_E8(99)` and
`D_E8(97)` dominates the corresponding negative-region bound.  Thus

```text
D_E8(99) >= 0,
D_E8(97) >= 0.
```

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound for
`n=99` and `n=97`.  Lemma 3 gives the negative-region bound in both cases.
Lemma 7 checks that the positive lower bound is larger than the negative
bound for both values of `n`.
