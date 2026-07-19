# DCT-RectTrig(E8, n=123) certificate

This note records the fixed wide-cap rectangular certificate for the single
E8 bridge step

```text
D_E8(123) = Q_3^{E8}(125) - 4 Q_3^{E8}(123).
```

## Lemma 1 (E8 constants and identities)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The root identities are

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800.
```

Proof.  These are Lemmas 1 and 3 of `DCT_RECT_TRIG_E8_N131.md`, replayed
by `direct_chain_rect_e8_n123_certificate.py`.

## Lemma 2 (wide-cap trigonometric bounds)

For `|z| <= 4`,

```text
2(1-cos z) >= z^2 - z^4/12 + z^6/480,
2 cos z >= 2 - z^2 + z^4/12 - z^6/360,
sin(z/2)/(z/2) >= exp(-z^2/20).
```

Proof.  Put `w=z^2`.  Since `w <= 16 < 90`, the Taylor tail for
`2cos z` after degree `6` starts with a positive term and is alternating
with decreasing terms.  Hence

```text
2 cos z >= 2 - z^2 + z^4/12 - z^6/360.
```

For the first bound, the alternating expansion gives

```text
2(1-cos z)
 >= z^2 - z^4/12 + z^6/360 - z^8/20160
    + z^10/1814400 - z^12/239500800.
```

It remains to check

```text
1/1440 - w/20160 + w^2/1814400 - w^3/239500800 >= 0
```

on `0 <= w <= 16`.  Its derivative is

```text
-1/20160 + w/907200 - w^2/79833600.
```

This derivative is increasing on the interval and at `w=16` equals
`-13/369600`, so the cubic is decreasing.  Its endpoint value is
`53/2138400 > 0`.

For the sine bound, put `x=z/2` and `y=x^2`.  Then `0 <= y <= 4`.  The
alternating Taylor expansion gives

```text
sin x / x >=
P(y) = 1 - y/6 + y^2/120 - y^3/5040
       + y^4/362880 - y^5/39916800.
```

Also

```text
exp(y/5) >= Q(y) = 1 + y/5 + y^2/50 + y^3/750.
```

A direct multiplication gives

```text
P(y)Q(y)-1 = y H(y)/29937600000,
```

where

```text
H(y) =
  997920000 - 149688000y - 15919200y^2 - 2768700y^3
  + 229590y^4 - 6420y^5 + 95y^6 - y^7.
```

On `[0,4]`, the Bernstein coefficients of `H` are

```text
997920000, 912384000, 5703033600/7, 4899037440/7,
3951015168/7, 2844846080/7, 1574170880/7, 19837696,
```

all positive.  Hence `H(y) >= 0`, so `P(y)Q(y) >= 1`.  Therefore

```text
sin x/x >= P(y) >= 1/Q(y) >= exp(-y/5),
```

which is the claimed sine bound.

## Lemma 3 (rational exponential lower bound)

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k) * (1 - r + r^2/2 - r^3/6).
```

Proof.  This is Lemma 3 of `DCT_RECT_TRIG_E8_N129.md`.

## Lemma 4 (wide-cap rectangular union)

Let `U_123` be the union of `1/40 x 1/40` cells

```text
R(s0,t0) = [s0,s0+1/40] x [t0,t0+1/40],
s0 in {53, 2121/40, ...},
t0 in {69/2, 1381/40, ..., 1499/40},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`28121` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 123 / (4*248) = 1845/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n123_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_123`,

```text
chi(y) - chi(x) >= (496/123) g_R,
g_R =
  s0 - t1
  - (248/(6*50*123)) s1^2
  + (1/480)(496^2/1800) s0^3/123^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  1 - (s1+t1)/123
  + (1/12)(496/50)(s0^2+t0^2)/123^2
  - (1/360)(496^2/1800)(s1^3+t1^3)/123^3.
```

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
Sine_R = ExpLower((4*248/(20*123))(s1+t1)),
Neg(123) = 520 * 16^123.
```

The rational lower ratio

```text
R123 =
  sum_R A0_lower * 496^123 * 123^(-250)
        * F_R * h_R^123
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(123)
```

satisfies `R123 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n123_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output
includes

```text
root_identity_constants: OK
e_upper_1457_536: OK
retained_cells=28121
log10_ratio123_lower=1.10
log10_best_cell_ratio123=-3.19
positive_gaps: OK
positive_h123: OK
positive_direct123: OK
ratio123_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=123))

The rectangular-union positive contribution to `D_E8(123)` dominates the
negative-region bound `520*16^123`.  Hence `D_E8(123) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=123`.  Lemma 7 checks that this lower bound is larger than the
negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=123)

Assume `E8-FiniteBridge(123)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 123`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(125) >= 4 Q_3^{E8}(123)`.  Theorem 8 of
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
