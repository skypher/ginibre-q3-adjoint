# DCT-RectTrig(E8, n=107) certificate

This note records the wide-cap rectangular certificate for the single E8
bridge step

```text
D_E8(107) = Q_3^{E8}(109) - 4 Q_3^{E8}(107).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities and cubic scalar trigonometric bounds are those of
Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker `direct_chain_rect_e8_n107_certificate.py` imports the
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

## Lemma 3 (asymmetric degree-28 negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character, set
`S = X+Y`, and set

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

On the negative-sign region for odd `n=107`,

```text
S in [-16,-2] union [0,2].
```

For every such `S`,

```text
|S^107(S^2-4)|
 <= (63/65) * 16^95 * S^12 * (S^2+4) * ((496-S)/512)^14.
```

Proof.  The right side is non-negative for all real `S`.  On
`S in [-16,-2]`, put `x=-S`.  After clearing the positive denominator
`65*512^14` and dividing by `x^12`, the residual is

```text
R_-(x) =
  63*16^95*(x^2+4)*(496+x)^14
  - 65*512^14*x^95*(x^2-4).
```

It has the endpoint factor `16-x`.  The checker divides by that factor and
verifies by exact Sturm root count that the quotient is positive on
`[2,16]`.  On `S in [0,2]`, after clearing the same denominator and
dividing by `S^12`, the residual is

```text
R_+(S) =
  63*16^95*(S^2+4)*(496-S)^14
  - 65*512^14*S^95*(4-S^2).
```

The checker verifies by exact Sturm root count that `R_+(S)` is positive on
`[0,2]`.  It returns

```text
asymmetric_negative_interval_sturm: OK
asymmetric_positive_interval_sturm: OK
```

The exact moment prefix `m_0,...,m_30` from
`character_ring_iter/logs/e8_70.log` and `references/oeis_A179663.txt`
gives

```text
J_12 = 5428655636,
J_14 = 1070740575644,
K_28 =
  sum_{i=0}^{14} binom(14,i) 496^{14-i}(-1)^i 512^{-14}
    * (J_{14+i}+4J_{12+i})
  = 9594616716049666093860048619540906615318349630451
    / 21267647932558653966460912964485513216.
```

Therefore the negative part of `D_E8(107)` is bounded by

```text
Neg(107) = (63/65) * 16^95 * K_28.
```

The checker returns

```text
moment_sources_0_30: OK
asymmetric_moment_integral: OK
```

## Lemma 4 (wide-cap rectangular union)

Let `U_107` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 517/10},
t0 in {30, 601/20, ..., 959/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`70562` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 107 / (4*248) = 1605/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n107_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_107`,

```text
chi(y) - chi(x) >= (496/107) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/107
       + L_3 (496^2/1800)s0^3/107^2
  - U_1 t1 - U_2 (496/50)t1^2/107
       - U_3 (496^2/1800)t1^3/107^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/107
  + C_2(496/50)(s0^2+t0^2)/107^2
  + C_3(496^2/1800)(s1^3+t1^3)/107^3.
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
    2(1/24) * 496(s1+t1)/107
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/107^2
  ),
Neg(107) = (63/65) * 16^95 * K_28.
```

The rational lower ratio

```text
R107 =
  sum_R A0_lower * 496^107 * 107^(-250)
        * F_R * h_R^107
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(107)
```

satisfies `R107 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n107_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=107 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=207/4 t_start=30 t_end=48
retained_cells=70562 cap=1605/31
asym_integral=9594616716049666093860048619540906615318349630451/21267647932558653966460912964485513216
log10_asym_negative_vs_p12=-0.40
negative_bound=(63/65)*16^95*asym_integral
log10_union_fraction_lower=-38.04
gap107_min=1428012189147/14311250000000000
h107_min=3400562485123005593/9494083250000000000 direct_factor107_min=11562359658241724095954587585909281649/366389246056884765625000000000000
sine_exponent_max=2032218299/50089375
log10_ratio107_lower=0.08
log10_best_cell_ratio107=-4.07
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
sine_quartic_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
moment_sources_0_30: OK
asymmetric_negative_interval_sturm: OK
asymmetric_positive_interval_sturm: OK
asymmetric_moment_integral: OK
cell_count_70562: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h107: OK
positive_direct107: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio107_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=107))

The rectangular-union positive contribution to `D_E8(107)` dominates the
negative-region bound `Neg(107)`.  Thus `D_E8(107) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=107`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that
the lower bound is larger than the negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=107)

Assume `E8-FiniteBridge(107)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 107`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(109) >= 4 Q_3^{E8}(107)`, hence the
assumed finite bridge through `n=107` gives `Q_3^{E8}(109) >= 0`.  This
supplies the hypothesis `E8-FiniteBridge(109)`, and Corollary 9 of
`DCT_RECT_TRIG_E8_N109.md` applies.
