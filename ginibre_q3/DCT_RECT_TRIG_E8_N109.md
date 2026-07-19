# DCT-RectTrig(E8, n=109) certificate

This note records the wide-cap rectangular certificate for the single E8
bridge step

```text
D_E8(109) = Q_3^{E8}(111) - 4 Q_3^{E8}(109).
```

## Lemma 1 (E8 constants, root identities, and scalar bounds)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities and cubic scalar trigonometric bounds are those of
Lemma 1 of `DCT_RECT_TRIG_E8_N111.md`.

Proof.  The checker `direct_chain_rect_e8_n109_certificate.py` imports the
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

## Lemma 3 (p=12 negative-region bound)

Let `X` and `Y` be independent copies of the E8 adjoint character, and set

```text
J_k = int int (X-Y)^2 (X+Y)^k dX dY.
```

The exact moment prefix `m_0,...,m_16` from
`character_ring_iter/logs/e8_70.log` and `references/oeis_A179663.txt`
gives

```text
J_12 = 5428655636,
J_14 = 1070740575644,
J_14 + 4J_12 = 1092455198188.
```

For odd `n >= 13`, the negative part of `D_E8(n)` is bounded by

```text
1092455198188 * 16^(n-12).
```

Proof.  The binomial identity for `J_k` is Lemma 5 of
`DCT_RECT_TRIG_E8_N113.md` with `k=12,14`.  On the negative-sign region,
with `S=X+Y`, one has `|S| <= 16`; therefore

```text
|S^n(S^2-4)| <= |S|^n(S^2+4)
              <= 16^(n-12) S^12(S^2+4).
```

Multiplying by `(X-Y)^2` and integrating gives the displayed bound.  The
checker returns

```text
moment_sources_0_16: OK
moment_negative_p12: OK
```

## Lemma 4 (wide-cap rectangular union)

Let `U_109` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 in {40, 801/20, ..., 1053/20},
t0 in {30, 601/20, ..., 959/20},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`77454` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 109 / (4*248) = 1635/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n109_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_109`,

```text
chi(y) - chi(x) >= (496/109) g_R,
g_R =
  L_1 s0 + L_2 (496/50)s0^2/109
       + L_3 (496^2/1800)s0^3/109^2
  - U_1 t1 - U_2 (496/50)t1^2/109
       - U_3 (496^2/1800)t1^3/109^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  (16 + 240 C_0)/496
  + C_1(s1+t1)/109
  + C_2(496/50)(s0^2+t0^2)/109^2
  + C_3(496^2/1800)(s1^3+t1^3)/109^3.
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
    2(1/24) * 496(s1+t1)/109
    + 2(1/2100) * (496^2/50)(s1^2+t1^2)/109^2
  ),
Neg(109) = 1092455198188 * 16^97.
```

The rational lower ratio

```text
R109 =
  sum_R A0_lower * 496^109 * 109^(-250)
        * F_R * h_R^109
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(109)
```

satisfies `R109 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n109_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=109 delta=4 e_upper=1457/536
grid_step=1/20 s_start=40 s_end=527/10 t_start=30 t_end=48
retained_cells=77454 cap=1635/31
negative_bound=1092455198188*16^97
log10_union_fraction_lower=-37.24
gap109_min=1262084876833/13366125000000000
h109_min=250022064633883337/690583125000000000 direct_factor109_min=62503278754418399015078220926255569/1938512342834472656250000000000
sine_exponent_max=52226353418/1299484375
log10_ratio109_lower=1.15
log10_best_cell_ratio109=-3.06
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
sine_quartic_sturm: OK
root_identity_constants: OK
e_upper_1457_536: OK
moment_sources_0_16: OK
moment_negative_p12: OK
cell_count_77454: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h109: OK
positive_direct109: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio109_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=109))

The rectangular-union positive contribution to `D_E8(109)` dominates the
negative-region bound `1092455198188*16^97`.  Thus `D_E8(109) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=109`.  Lemma 3 gives the negative-region bound.  Lemma 7 checks that
the lower bound is larger than the negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=109)

Assume `E8-FiniteBridge(109)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 109`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(111) >= 4 Q_3^{E8}(109)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N111.md` gives
`Q_3^{E8}(113) >= 4 Q_3^{E8}(111)`.  Corollary 9 of
`DCT_RECT_TRIG_E8_N111.md` then applies with the assumed finite bridge.
