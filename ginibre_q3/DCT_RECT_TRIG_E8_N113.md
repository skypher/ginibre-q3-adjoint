# DCT-RectTrig(E8, n=113) certificate

This note records the wide-cap rectangular certificate for the single E8
bridge step

```text
D_E8(113) = Q_3^{E8}(115) - 4 Q_3^{E8}(113).
```

## Lemma 1 (E8 constants and identities)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

With `Q(u)=sum_{alpha>0}(alpha.u)^2`,

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800.
```

Proof.  These are Lemmas 1 and 3 of `DCT_RECT_TRIG_E8_N131.md`, replayed
by `direct_chain_rect_e8_n113_certificate.py` through the shared
`root_identity_check()`.

## Lemma 2 (cubic scalar trigonometric bounds)

For `0 <= w <= 16`, set

```text
L(w) =
  (99119225/100000000)w
  - (7802858/100000000)w^2
  + (181097/100000000)w^3,

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
L(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  This is Lemma 2 of `DCT_RECT_TRIG_E8_N117.md`, replayed in
`direct_chain_rect_e8_n113_certificate.py` by exact Sturm root counts.  The
checker returns

```text
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
```

## Lemma 3 (quartic Weyl-sine bound)

For `|z| <= 4`,

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100.
```

Proof.  This is Lemma 3 of `DCT_RECT_TRIG_E8_N115.md`, replayed in
`direct_chain_rect_e8_n113_certificate.py` by the same exact Sturm root
count after the substitutions `x=z/2` and `y=x^2`.  The checker returns

```text
sine_quartic_sturm: OK
```

## Lemma 4 (rational exponential lower bound)

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k)
  * (1 - r + r^2/2 - r^3/6 + r^4/24
     - r^5/120 + r^6/720 - r^7/5040).
```

Proof.  This is Lemma 3 of `DCT_RECT_TRIG_E8_N121.md`, replayed by the
checker through `e_upper_1457_536: OK`.

## Lemma 5 (low-moment negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character.  The
exact moment prefix is

```text
m_0 = 1, m_1 = 0, m_2 = 1, m_3 = 1,
m_4 = 5, m_5 = 16, m_6 = 79.
```

For

```text
J_k = int int (X-Y)^2 (X+Y)^k dX dY,
```

the binomial identity

```text
J_k =
  sum_{i=0}^k binom(k,i)
    (m_{i+2}m_{k-i} + m_i m_{k-i+2}
     - 2m_{i+1}m_{k-i+1})
```

gives

```text
J_2 = 8,   J_4 = 144,   J_4 + 4J_2 = 176.
```

For odd `n >= 3`, the negative part of `D_E8(n)` is bounded by

```text
176 * 16^(n-2).
```

Proof.  The moment values are checked against both
`character_ring_iter/logs/e8_70.log` and `references/oeis_A179663.txt` by
`low_moment_sources: OK`.  The displayed binomial identity follows by
expanding `(X-Y)^2(X+Y)^k` and using independence.  On the negative-sign
region for `D_E8(n)`, writing `S=X+Y`, either `S<0` or `0<S<2`; since
`min chi_adj = -8`, this gives `|S| <= 16`.  Therefore

```text
|S^n(S^2-4)| <= |S|^n(S^2+4)
              <= 16^(n-2) S^2(S^2+4).
```

Multiplying by `(X-Y)^2` and integrating gives
`16^(n-2)(J_4+4J_2) = 176*16^(n-2)`.  The checker returns

```text
low_moment_negative_176: OK
```

## Lemma 6 (wide-cap rectangular union)

Let `U_113` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 273/5},
t0 in {30, 601/20, ..., 959/20},
```

retained by the positive rational gap in Lemma 7.  The exact enumeration has
`91581` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 113 / (4*248) = 1695/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n113_certificate.py`.

## Lemma 7 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_113`,

```text
chi(y) - chi(x) >= (496/113) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/113
       + L_3 (496^2/1800)s0^3/113^2
  - U_1 t1 - U_2 (496/50)t1^2/113
       - U_3 (496^2/1800)t1^3/113^2,
```

where `(L_1,L_2,L_3)` and `(U_1,U_2,U_3)` are the coefficients of `L` and
`U` from Lemma 2.  Also

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/113
  + C_2(496/50)(s0^2+t0^2)/113^2
  + C_3(496^2/1800)(s1^3+t1^3)/113^3,
```

where `(C_0,C_1,C_2,C_3)` are the coefficients of `C` from Lemma 2.
Moreover every retained cell has `g_R > 0`, `h_R > 0`, and
`496^2 h_R^2 - 4 > 0`.

Proof.  Lemmas 1 and 2 give the displayed formulas.  The endpoint choices
use lower endpoints in positive terms and upper endpoints in negative terms.
The exact rational checker verifies the positivity comparisons over all
retained cells.

## Proposition 8 (rectangular radial mass lower)

For a retained cell `R=[s0,s1] x [t0,t1]`, define

```text
F_R =
  g_R^2 * ExpLower(s1+t1)
  * (s1^124-s0^124)/124
  * (t1^124-t0^124)/124
  * 1/(123!)^2
  / 124,
```

where `ExpLower` is the rational lower bound from Lemma 4.  This is a lower
bound for the normalized radial mass factor on `R`.

Proof.  On the cell the gap factor is at least `g_R^2`; Lemma 4 supplies
the radial Gaussian lower bound.  Since `a=124`, `Gamma(124)=123!`, giving
the displayed factor.

## Lemma 9 (exact-rational arithmetic check)

Let

```text
A0_lower = A0_prefactor(E8) * (113/355)^8,
Sine_R =
  ExpLower(
    2(1/24) * 496(s1+t1)/113
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/113^2
  ),
Neg(113) = 176 * 16^111.
```

The rational lower ratio

```text
R113 =
  sum_R A0_lower * 496^113 * 113^(-250)
        * F_R * h_R^113
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(113)
```

satisfies `R113 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n113_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=113 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=1093/20 t_start=30 t_end=48
retained_cells=91581 cap=1695/31
negative_bound=176*16^111
log10_union_fraction_lower=-35.78
gap113_min=30644144709323/143651250000000000
h113_min=61950654956354875633/167736776250000000000 direct_factor113_min=3837426189419445291744113611010897150689/114365025472906494140625000000000000
sine_exponent_max=165452401781/4189828125
log10_ratio113_lower=1.89
log10_best_cell_ratio113=-2.44
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
sine_quartic_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
low_moment_sources: OK
low_moment_negative_176: OK
cell_count_91581: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h113: OK
positive_direct113: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio113_gt_1: OK
```

## Theorem 10 (DCT-RectTrig(E8, n=113))

The rectangular-union positive contribution to `D_E8(113)` dominates the
negative-region bound `176*16^111`.  Thus `D_E8(113) >= 0`.

Proof.  Lemmas 6-8 give the positive rectangular-union lower bound at
`n=113`.  Lemma 5 gives the negative-region bound.  Lemma 9 checks that
the lower bound is larger than the negative-region bound.

## Corollary 11 (conditional E8 odd tail from n=113)

Assume `E8-FiniteBridge(113)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 113`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 10 gives `Q_3^{E8}(115) >= 4 Q_3^{E8}(113)`.  Theorem 9 of
`DCT_RECT_TRIG_E8_N115.md` gives
`Q_3^{E8}(117) >= 4 Q_3^{E8}(115)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N117.md` gives
`Q_3^{E8}(119) >= 4 Q_3^{E8}(117)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N119.md` gives
`Q_3^{E8}(121) >= 4 Q_3^{E8}(119)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N121.md` gives
`Q_3^{E8}(123) >= 4 Q_3^{E8}(121)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N123.md` gives
`Q_3^{E8}(125) >= 4 Q_3^{E8}(123)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N125.md` gives
`Q_3^{E8}(127) >= 4 Q_3^{E8}(125)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N127.md` gives
`Q_3^{E8}(129) >= 4 Q_3^{E8}(127)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N129.md` gives
`Q_3^{E8}(131) >= 4 Q_3^{E8}(129)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N131.md` gives
`Q_3^{E8}(133) >= 4 Q_3^{E8}(131)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_TAIL.md` gives
`Q_3^{E8}(n+2) >= 4 Q_3^{E8}(n)` for every odd `n >= 133`.  Induction over
odd indices from the assumed finite bridge gives the claim.
