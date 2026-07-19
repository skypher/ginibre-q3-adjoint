# DCT-RectTrig(E8, n=93) certificate

This note records the wide-cap rectangular certificate for the E8 bridge step

```text
D_E8(93) = Q_3^{E8}(95) - 4 Q_3^{E8}(93).
```

## Lemma 1 (E8 constants, root identities, and widened scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800
```

are those of Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

For `0 <= w <= (21/5)^2`, set

```text
L_*(w) =
  (99119225/100000000)w
  - (7802858/100000000)w^2
  + (7/4000)w^3,

U(w) =
  w
  - (8228444/100000000)w^2
  + (233397/100000000)w^3,

C(w) =
  199909270/100000000
  - (99730681/100000000)w
  + (8092566/100000000)w^2
  - (216849/100000000)w^3.
```

Then

```text
L_*(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `direct_chain_rect_e8_n93_certificate.py` verifies the
root identities and checks the three residual polynomials by exact Sturm
root counts on `[0,(21/5)^2]`.  It returns

```text
wide_gap_lower_sturm: OK
wide_gap_upper_sturm: OK
wide_cos_lower_sturm: OK
root_identity_constants: OK
```

## Lemma 2 (widened quartic Weyl-sine and exponential bounds)

For `|z| <= 21/5`,

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2000.
```

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k)
  * (1 - r + r^2/2 - r^3/6 + r^4/24
     - r^5/120 + r^6/720 - r^7/5040).
```

Proof.  With `y=(z/2)^2`, the sine bound is checked by exact Sturm
positivity for the residual of

```text
sin(sqrt(y))/sqrt(y) * exp(y/6 + y^2/125) - 1
```

after replacing both factors by the same Taylor-side polynomials used in the
earlier E8 bridge checks.  The exponential bound is the same rational
fractional-part bound as in `DCT_RECT_TRIG_E8_N111.md`.  The checker returns

```text
sine_quartic_21_5_sturm: OK
e_upper_1457_536: OK
```

## Lemma 3 (shifted Chebyshev degree-27 negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character, set
`S = X+Y`, and set

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

Let

```text
C_27(S) = T_27((2S-503)/489) / T_27(-535/489),
```

where `T_27` is the Chebyshev polynomial of the first kind.  Then
`C_27(-16)=1`.  On the negative-sign region for odd `n=93`,

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^93(S^2-4)|
 <= (63/65) * 16^81 * S^12 * (S^2+4) * C_27(S)^2.
```

Proof.  The right side is non-negative for all real `S`.  On
`S in [-16,-2]`, put `x=-S`.  After clearing positive denominators and
dividing by `x^12(16-x)`, the checker verifies positivity on `[2,16]` by
checking that every Bernstein coefficient of the residual quotient is
positive.  On `S in [0,2]`, after clearing positive denominators and
dividing by `S^12`, the checker verifies positivity on `[0,2]` by the same
Bernstein-coefficient test.  It returns

```text
shifted_chebyshev_normalization: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
```

The exact moment prefix `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log`, with OEIS overlap `m_0,...,m_30`
checked against `references/oeis_A179663.txt`, gives the same integral as
the `n=95` certificate:

```text
K_70 =
  int int (X-Y)^2 S^12(S^2+4) C_27(S)^2 dX dY
  = 13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828
    / 51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025.
```

Therefore the negative part of `D_E8(93)` is bounded by

```text
Neg(93) = (63/65) * 16^81 * K_70.
```

The checker returns

```text
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_moment_integral: OK
```

## Lemma 4 (widened rectangular union)

Let `U_93` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 991/20},
t0 in {30, 601/20, ..., 242/5},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`53822` cells.  Each cell lies inside the root cap `|alpha.v| <= 21/5`
guaranteed by

```text
s,t <= (21/5)^2 * 30 * 93 / (4*248) = 3969/80.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n93_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_93`,

```text
chi(y) - chi(x) >= (496/93) g_R,
g_R =
  L_{*,1} s0 + L_{*,2} (496/50)s0^2/93
       + L_{*,3} (496^2/1800)s0^3/93^2
  - U_1 t1 - U_2 (496/50)t1^2/93
       - U_3 (496^2/1800)t1^3/93^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/93
  + C_2(496/50)(s0^2+t0^2)/93^2
  + C_3(496^2/1800)(s1^3+t1^3)/93^3.
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
    2(1/24) * 496(s1+t1)/93
    + 2(1/2000) * (496^2/50)(s1^2+t1^2)/93^2
  ),
Neg(93) = (63/65) * 16^81 * K_70.
```

The rational lower ratio

```text
R93 =
  sum_R A0_lower * 496^93 * 93^(-250)
        * F_R * h_R^93
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(93)
```

satisfies `R93 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n93_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=93 delta=21/5 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=248/5 t_start=30 t_end=969/20
retained_cells=53822 cap=3969/80
wide_gap_lower_coeffs=99119225/100000000,-7802858/100000000,7/4000
sine_bound=exp(-z^2/24-z^4/2000)
shifted_chebyshev_degree=27 shifted_interval=[7,496] moment_max=70 cheb_at_0=1286950147423114883724792415855906842289505985166769366698842399977044845479/227714441580758816060790539068573206609030282743660865709718401498816010950855
cheb_integral=13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828/51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025
log10_cheb_negative_vs_p12=-9.65
negative_bound=(63/65)*16^81*cheb_integral
log10_union_fraction_lower=-40.58
gap93_min=7225097/40500000000
h93_min=962182419750323/3138750000000000 direct_factor93_min=925634828212524259905658604329/40045166015625000000000000
sine_exponent_max=434182/9375
log10_ratio93_lower=0.04
log10_best_cell_ratio93=-4.04
wide_gap_lower_sturm: OK
wide_gap_upper_sturm: OK
wide_cos_lower_sturm: OK
sine_quartic_21_5_sturm: OK
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
shifted_chebyshev_moment_integral: OK
ratio93_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=93))

The rectangular-union positive contribution to `D_E8(93)` dominates the
negative-region bound `Neg(93)`.  Thus

```text
D_E8(93) >= 0.
```

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=93`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that the
positive lower bound is larger than the negative bound.
