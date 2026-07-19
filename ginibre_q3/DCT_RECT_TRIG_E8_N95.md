# DCT-RectTrig(E8, n=95) certificate

This note records the wide-cap rectangular certificate for the E8 bridge step

```text
D_E8(95) = Q_3^{E8}(97) - 4 Q_3^{E8}(95).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities and cubic scalar trigonometric bounds are those of
Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker `direct_chain_rect_e8_n95_certificate.py` imports the
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
`C_27(-16)=1`.  On the negative-sign region for odd `n=95`,

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^95(S^2-4)|
 <= (63/65) * 16^83 * S^12 * (S^2+4) * C_27(S)^2.
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
checked against `references/oeis_A179663.txt`, gives

```text
K_70 =
  int int (X-Y)^2 S^12(S^2+4) C_27(S)^2 dX dY
  = 13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828
    / 51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025.
```

Therefore the negative part of `D_E8(95)` is bounded by

```text
Neg(95) = (63/65) * 16^83 * K_70.
```

The checker returns

```text
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_moment_integral: OK
```

## Lemma 4 (wide-cap rectangular union)

Let `U_95` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 459/10},
t0 in {30, 601/20, ..., 453/10},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`29784` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 15*95/31 = 1425/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n95_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_95`,

```text
chi(y) - chi(x) >= (496/95) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/95
       + L_3 (496^2/1800)s0^3/95^2
  - U_1 t1 - U_2 (496/50)t1^2/95
       - U_3 (496^2/1800)t1^3/95^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/95
  + C_2(496/50)(s0^2+t0^2)/95^2
  + C_3(496^2/1800)(s1^3+t1^3)/95^3.
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
    2(1/24) * 496(s1+t1)/95
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/95^2
  ),
Neg(95) = (63/65) * 16^83 * K_70.
```

The rational lower ratio

```text
R95 =
  sum_R A0_lower * 496^95 * 95^(-250)
        * F_R * h_R^95
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(95)
```

satisfies `R95 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n95_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=95 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=919/20 t_start=30 t_end=907/20
retained_cells=29784 cap=1425/31
shifted_chebyshev_degree=27 shifted_interval=[7,496] moment_max=70 cheb_at_0=1286950147423114883724792415855906842289505985166769366698842399977044845479/227714441580758816060790539068573206609030282743660865709718401498816010950855
cheb_integral=13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828/51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025
log10_cheb_negative_vs_p12=-9.65
negative_bound=(63/65)*16^83*cheb_integral
log10_union_fraction_lower=-45.39
gap95_min=9482346092687/101531250000000000
h95_min=17205807832136089069/49834921875000000000 direct_factor95_min=295999443352833420655877472050701286761/10094950890541076660156250000000000
sine_exponent_max=24808619798/592265625
log10_ratio95_lower=0.70
log10_best_cell_ratio95=-3.01
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
shifted_chebyshev_moment_integral: OK
ratio95_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=95))

The rectangular-union positive contribution to `D_E8(95)` dominates the
negative-region bound `Neg(95)`.  Thus

```text
D_E8(95) >= 0.
```

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=95`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that
the positive lower bound is larger than the negative bound.
