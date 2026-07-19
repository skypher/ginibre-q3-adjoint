# DCT-RectTrig(E8, n=117) certificate

This note records the fixed wide-cap rectangular certificate for the single
E8 bridge step

```text
D_E8(117) = Q_3^{E8}(119) - 4 Q_3^{E8}(117).
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
by `direct_chain_rect_e8_n117_certificate.py` through the shared
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

Proof.  The alternating Taylor expansions give, on `0 <= w <= 16`,

```text
2(1-cos sqrt(w))
 >= w - w^2/12 + w^3/360 - w^4/20160
    + w^5/1814400 - w^6/239500800,

2(1-cos sqrt(w))
 <= w - w^2/12 + w^3/360 - w^4/20160 + w^5/1814400,

2 cos sqrt(w)
 >= 2 - w + w^2/12 - w^3/360 + w^4/20160
    - w^5/1814400 + w^6/239500800 - w^7/43589145600.
```

Subtracting the displayed cubic polynomials leaves rational polynomials.
After removing the endpoint factors `w` and `w^2` in the first two
residuals, `direct_chain_rect_e8_n117_certificate.py` checks the three
remaining residuals by exact Sturm root counts on `[0,16]`.  The checker
returns

```text
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
```

## Lemma 3 (wide-cap sine and exponential bounds)

For `|z| <= 4`,

```text
sin(z/2)/(z/2) >= exp(-4z^2/81).
```

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k)
  * (1 - r + r^2/2 - r^3/6 + r^4/24
     - r^5/120 + r^6/720 - r^7/5040).
```

Proof.  The sine bound is Lemma 2 of `DCT_RECT_TRIG_E8_N121.md`.  The
exponential bound is Lemma 3 of `DCT_RECT_TRIG_E8_N121.md`, and the checker
again verifies `e < 1457/536`.

## Lemma 4 (wide-cap rectangular union)

Let `U_117` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 1131/20},
t0 in {30, 601/20, ..., 999/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`111427` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 117 / (4*248) = 1755/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n117_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_117`,

```text
chi(y) - chi(x) >= (496/117) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/117
       + L_3 (496^2/1800)s0^3/117^2
  - U_1 t1 - U_2 (496/50)t1^2/117
       - U_3 (496^2/1800)t1^3/117^2,
```

where `(L_1,L_2,L_3)` and `(U_1,U_2,U_3)` are the coefficients of `L` and
`U` from Lemma 2.  Also

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/117
  + C_2(496/50)(s0^2+t0^2)/117^2
  + C_3(496^2/1800)(s1^3+t1^3)/117^3,
```

where `(C_0,C_1,C_2,C_3)` are the coefficients of `C` from Lemma 2.
Moreover every retained cell has `g_R > 0`, `h_R > 0`, and
`496^2 h_R^2 - 4 > 0`.

Proof.  Lemmas 1 and 2 give the displayed formulas.  The endpoint choices
use lower endpoints in positive terms and upper endpoints in negative terms.
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
  / 124,
```

where `ExpLower` is the rational lower bound from Lemma 3.  This is a lower
bound for the normalized radial mass factor on `R`.

Proof.  On the cell the gap factor is at least `g_R^2`; Lemma 3 supplies
the radial Gaussian lower bound.  Since `a=124`, `Gamma(124)=123!`, giving
the displayed factor.

## Lemma 7 (exact-rational arithmetic check)

Let

```text
A0_lower = A0_prefactor(E8) * (113/355)^8,
Sine_R = ExpLower((4*248/117)/(81/4) * (s1+t1)),
Neg(117) = 520 * 16^117.
```

The rational lower ratio

```text
R117 =
  sum_R A0_lower * 496^117 * 117^(-250)
        * F_R * h_R^117
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(117)
```

satisfies `R117 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n117_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=117 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=283/5 t_start=30 t_end=50
retained_cells=111427 cap=1755/31
log10_union_fraction_lower=-33.40
gap117_min=183481292003/3422250000000000
h117_min=1317904655485969331/3580529062500000000 direct_factor117_min=1736664236161538716850281550672587561/52111197513198852539062500000000
sine_exponent_max=162688/3645
log10_ratio117_lower=0.05
log10_best_cell_ratio117=-4.38
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
cell_count_111427: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h117: OK
positive_direct117: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio117_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=117))

The rectangular-union positive contribution to `D_E8(117)` dominates the
negative-region bound `520*16^117`.  Hence `D_E8(117) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=117`.  Lemma 7 checks that this lower bound is larger than the
negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=117)

Assume `E8-FiniteBridge(117)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 117`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(119) >= 4 Q_3^{E8}(117)`.  Theorem 8 of
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
