# DCT-RectTrig(E8, n=111) certificate

This note records the wide-cap rectangular certificate for the single E8
bridge step

```text
D_E8(111) = Q_3^{E8}(113) - 4 Q_3^{E8}(111).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800
```

and the cubic scalar bounds

```text
L(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w),             0 <= w <= 16,
```

are those of Lemmas 1 and 2 of `DCT_RECT_TRIG_E8_N113.md`.

Proof.  The checker `direct_chain_rect_e8_n111_certificate.py` imports the
same exact Sturm and root-identity checks and returns

```text
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
root_identity_constants: OK
```

## Lemma 2 (quartic Weyl-sine and exponential bounds)

For `|z| <= 4`,

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100.
```

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k)
  * (1 - r + r^2/2 - r^3/6 + r^4/24
     - r^5/120 + r^6/720 - r^7/5040).
```

Proof.  These are Lemmas 3 and 4 of `DCT_RECT_TRIG_E8_N113.md`, replayed by
the checker.  It returns

```text
sine_quartic_sturm: OK
e_upper_1457_536: OK
```

## Lemma 3 (low-moment negative-region bound)

For odd `n >= 3`, the negative part of `D_E8(n)` is bounded by

```text
176 * 16^(n-2).
```

Proof.  This is Lemma 5 of `DCT_RECT_TRIG_E8_N113.md`.  The checker verifies
the E8 moment prefix against both local exact moments and OEIS A179663, and
then checks the binomial expansion

```text
int int (X-Y)^2 (X+Y)^2 ((X+Y)^2+4) dX dY = 176.
```

It returns

```text
low_moment_sources: OK
low_moment_negative_176: OK
```

## Lemma 4 (wide-cap rectangular union)

Let `U_111` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 1073/20},
t0 in {30, 601/20, ..., 959/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`84702` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 111 / (4*248) = 1665/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n111_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_111`,

```text
chi(y) - chi(x) >= (496/111) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/111
       + L_3 (496^2/1800)s0^3/111^2
  - U_1 t1 - U_2 (496/50)t1^2/111
       - U_3 (496^2/1800)t1^3/111^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/111
  + C_2(496/50)(s0^2+t0^2)/111^2
  + C_3(496^2/1800)(s1^3+t1^3)/111^3.
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
    2(1/24) * 496(s1+t1)/111
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/111^2
  ),
Neg(111) = 176 * 16^109.
```

The rational lower ratio

```text
R111 =
  sum_R A0_lower * 496^111 * 111^(-250)
        * F_R * h_R^111
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(111)
```

satisfies `R111 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n111_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=111 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=537/10 t_start=30 t_end=48
retained_cells=84702 cap=1665/31
negative_bound=176*16^109
log10_union_fraction_lower=-36.47
gap111_min=11131456399439/138611250000000000
h111_min=358832265265913163/981401875000000000 direct_factor111_min=128744934644639863575293200256664569/3914987806701660156250000000000
sine_exponent_max=17897731766/449203125
log10_ratio111_lower=0.39
log10_best_cell_ratio111=-3.88
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
sine_quartic_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
low_moment_sources: OK
low_moment_negative_176: OK
cell_count_84702: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h111: OK
positive_direct111: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio111_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=111))

The rectangular-union positive contribution to `D_E8(111)` dominates the
negative-region bound `176*16^109`.  Thus `D_E8(111) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=111`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that
the lower bound is larger than the negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=111)

Assume `E8-FiniteBridge(111)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 111`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(113) >= 4 Q_3^{E8}(111)`.  Theorem 10 of
`DCT_RECT_TRIG_E8_N113.md` gives
`Q_3^{E8}(115) >= 4 Q_3^{E8}(113)`.  Corollary 11 of
`DCT_RECT_TRIG_E8_N113.md` then applies with the assumed finite bridge.
