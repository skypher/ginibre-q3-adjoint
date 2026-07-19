# DCT-RectTrig(E8, n=121) certificate

This note records the fixed wide-cap rectangular certificate for the single
E8 bridge step

```text
D_E8(121) = Q_3^{E8}(123) - 4 Q_3^{E8}(121).
```

## Lemma 1 (E8 constants and root identities)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

With `Q(u)=sum_{alpha>0}(alpha.u)^2`,

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800,
sum_{alpha>0}(alpha.u)^8 >= Q(u)^4/300000.
```

Proof.  The constants and the fourth- and sixth-power identities are Lemmas
1 and 3 of `DCT_RECT_TRIG_E8_N131.md`, replayed by
`direct_chain_rect_e8_n121_certificate.py` through the shared
`root_identity_check()`.

For the eighth-power lower bound, put `a_alpha=(alpha.u)^4`.  Cauchy's
inequality over the `120` positive roots gives

```text
sum_{alpha>0}(alpha.u)^8
  >= (sum_{alpha>0}(alpha.u)^4)^2/120
   = (Q(u)^2/50)^2/120
   = Q(u)^4/300000.
```

## Lemma 2 (wide-cap trigonometric bounds)

For `|z| <= 4`,

```text
2(1-cos z) >= z^2 - z^4/12 + z^6/475,
2 cos z >= 2 - z^2 + z^4/12 - z^6/360 + z^8/30000,
sin(z/2)/(z/2) >= exp(-4z^2/81).
```

Proof.  Put `w=z^2`.  Since `w <= 16 < 90`, the Taylor tails used below
are alternating with decreasing terms.

For the first bound,

```text
2(1-cos z)
 >= z^2 - z^4/12 + z^6/360 - z^8/20160
    + z^10/1814400 - z^12/239500800.
```

It remains to check

```text
1/360 - 1/475 - w/20160 + w^2/1814400 - w^3/239500800 >= 0
```

on `0 <= w <= 16`.  Its derivative is

```text
-1/20160 + w/907200 - w^2/79833600.
```

This derivative is increasing on the interval and at `w=16` equals
`-13/369600`, so the cubic is decreasing.  Its endpoint value is
`29/10157400 > 0`.

For the second bound, the alternating expansion gives

```text
2 cos z
 >= 2 - z^2 + z^4/12 - z^6/360
    + z^8/20160 - z^10/1814400.
```

It remains to check

```text
1/20160 - 1/30000 - w/1814400 >= 0
```

on `0 <= w <= 16`.  The left side is decreasing and its endpoint value is
`169/22680000 > 0`.

For the sine bound, put `x=z/2` and `y=x^2`, so `0 <= y <= 4`.  The
alternating Taylor expansion gives

```text
sin x / x >=
P(y) = 1 - y/6 + y^2/120 - y^3/5040
       + y^4/362880 - y^5/39916800.
```

Also

```text
exp(16y/81) >= Q(y) = sum_{j=0}^4 (16y/81)^j/j!.
```

A direct multiplication gives

```text
P(y)Q(y)-1 = y H(y)/5154862058438400,
```

where

```text
H(y) =
  159100680816000 - 26182865126880y - 2676926316240y^2
  - 126387254070y^3 - 16595723043y^4 + 1662781392y^5
  - 49152384y^6 + 735232y^7 - 8192y^8.
```

On `[0,4]`, the degree-eight Bernstein coefficients of `H` are

```text
159100680816000,
146009248252560,
919717004558880/7,
805650468052320/7,
3392145898528896/35,
536405176485072/7,
53982997834464,
28747476133184,
713979039872,
```

all positive.  Hence `H(y) >= 0`, so `P(y)Q(y) >= 1`.  Therefore

```text
sin x/x >= P(y) >= 1/Q(y) >= exp(-16y/81) = exp(-4z^2/81).
```

## Lemma 3 (rational exponential lower bound)

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k)
  * (1 - r + r^2/2 - r^3/6 + r^4/24
     - r^5/120 + r^6/720 - r^7/5040).
```

Proof.  The checker verifies

```text
sum_{j=0}^7 1/j! + 1/(7!*7) < 1457/536,
```

so `e < 1457/536`.  The alternating Taylor expansion for `exp(-r)` on
`0 <= r < 1` gives the displayed degree-seven lower polynomial.

## Lemma 4 (wide-cap rectangular union)

Let `U_121` be the union of `1/30 x 1/30` cells

```text
R(s0,t0) = [s0,s0+1/30] x [t0,t0+1/30],
s0 in {45, 1351/30, ..., 117/2},
t0 in {30, 901/30, ..., 1259/30},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`74053` cells.  Each cell lies inside the root cap `|alpha.v| <= 4`
guaranteed by

```text
s,t <= 4^2 * 30 * 121 / (4*248) = 1815/31.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n121_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_121`,

```text
chi(y) - chi(x) >= (496/121) g_R,
g_R =
  s0 - t1
  - (248/(6*50*121)) s1^2
  + (1/475)(496^2/1800) s0^3/121^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  1 - (s1+t1)/121
  + (1/12)(496/50)(s0^2+t0^2)/121^2
  - (1/360)(496^2/1800)(s1^3+t1^3)/121^3
  + (1/30000)(496^3/300000)(s0^4+t0^4)/121^4.
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
Sine_R = ExpLower((4*248/121)/(81/4) * (s1+t1)),
Neg(121) = 520 * 16^121.
```

The rational lower ratio

```text
R121 =
  sum_R A0_lower * 496^121 * 121^(-250)
        * F_R * h_R^121
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(121)
```

satisfies `R121 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n121_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output is

```text
group=E8 n=121 delta=4 e_upper=1457/536
grid_step=1/30 s_start=45 s_end=60 t_start=30 t_end=42
retained_cells=74053 cap=1815/31
log10_union_fraction_lower=-42.58
gap121_min=89528719/1173567656250
h121_min=239264176093148353904/572269131771240234375 direct_factor121_min=14082453096235815349484722043995321859725232956/327491959178209119030846655368804931640625
sine_exponent_max=527744/13365
log10_ratio121_lower=0.10
log10_best_cell_ratio121=-4.29
root_identity_constants: OK
e_upper_1457_536: OK
cell_count_74053: OK
all_cells_in_cap: OK
positive_gaps: OK
positive_h121: OK
positive_direct121: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
ratio121_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=121))

The rectangular-union positive contribution to `D_E8(121)` dominates the
negative-region bound `520*16^121`.  Hence `D_E8(121) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=121`.  Lemma 7 checks that this lower bound is larger than the
negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=121)

Assume `E8-FiniteBridge(121)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 121`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(123) >= 4 Q_3^{E8}(121)`.  Theorem 8 of
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
