# DCT-RectTrig(E8, n=131) certificate

This note records the fixed wide-cap rectangular certificate for the single
E8 bridge step

```text
D_E8(131) = Q_3^{E8}(133) - 4 Q_3^{E8}(131).
```

## Lemma 1 (E8 constants)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

Proof.  These are the E8 entries in the direct-tail root-data table.  They
are replayed by `direct_chain_rect_e8_n131_certificate.py` from the same
`ROOT_DATA` table used by the A0 and truncated-mass checkers.

## Lemma 2 (wide-cap trigonometric bounds)

For `|z| <= 19/6`,

```text
2(1-cos z) >= z^2 - z^4/12 + z^6/432,
2 cos z >= 2 - z^2 + z^4/12 - z^6/360,
sin(z/2)/(z/2) >= exp(-z^2/21).
```

Proof.  Put `w=z^2`.  Since `w <= 361/36 < 11`, the Taylor tails for
`cos z` have decreasing alternating terms from degree `8` onward.  Thus

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
1/2160 - w/20160 + w^2/1814400 - w^3/239500800 >= 0
```

on `0 <= w <= 361/36`.  Its derivative is

```text
-1/20160 + w/907200 - w^2/79833600.
```

The derivative is increasing on this interval and at `w=361/36` equals
`-4118833/103464345600`, so the cubic is decreasing.  Its endpoint value is
`26761073/1596307046400 > 0`.

For the sine bound, put `x=z/2` and

```text
F(x) = log(sin x/x) + 4x^2/21.
```

On `0 <= x <= pi/2`, the partial fraction identity

```text
cot x = 1/x - 2x sum_{m>=1} 1/(m^2 pi^2 - x^2)
```

gives

```text
F'(x)/x = 8/21 - 2 sum_{m>=1} 1/(m^2 pi^2 - x^2).
```

The sum is increasing, so the minimum of `F` on this interval is at an
endpoint.  At `0`, `F=0`.  At `pi/2`, the inequality is
`log(pi/2) < pi^2/21`, proved in Lemma 2 of
`DCT_RECT_TRIG_E8_TAIL.md`.

On `pi/2 <= x <= 19/12`, `cot x <= 0` and
`8x^2/21 <= 361/378 < 1`, hence

```text
F'(x) <= -1/x + 8x/21 < 0.
```

So the minimum on this interval is at `x=19/12`.  Using
`pi < 22/7 < 19/6` and `pi > 333/106`, write
`19/12 = pi/2 + y` with `0 <= y < 2/159`; then

```text
sin(19/12)/(19/12) = cos y /(19/12)
  > (12/19)(1 - (2/159)^2/2)
  = 101116/160113.
```

For `A=361/756`,

```text
exp(-A) <= 1/(1 + A + A^2/2),
```

and the exact rational comparison

```text
(101116/160113)(1 + A + A^2/2) - 1
 = 4396547/863305128 > 0
```

gives `F(19/12) > 0`.

## Lemma 3 (E8 quartic and sextic root identities)

For E8, with `Q(u)=sum_{alpha>0}(alpha.u)^2`,

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50,
sum_{alpha>0}(alpha.u)^6 = Q(u)^3/1800.
```

Proof.  The E8 invariant degrees are `2,8,12,14,18,20,24,30`, so the
degree `4` and degree `6` Weyl-invariant polynomials above are scalar
multiples of `Q^2` and `Q^3`.  Evaluating at a root `beta`, the exact root
enumeration in `direct_chain_rect_e8_n131_certificate.py` gives

```text
sum_{alpha>0}(alpha.beta)^2 = 60,
sum_{alpha>0}(alpha.beta)^4 = 72,
sum_{alpha>0}(alpha.beta)^6 = 120.
```

Thus the scalars are `72/60^2 = 1/50` and `120/60^3 = 1/1800`.

## Lemma 4 (wide-cap rectangular union)

Let `U_131` be the union of `1/40 x 1/40` cells

```text
R(s0,t0) = [s0,s0+1/40] x [t0,t0+1/40],
s0 in {37, 1481/40, ...},
t0 in {26, 1041/40, ..., 1239/40},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`17797` cells.  Each cell lies inside the root cap `|alpha.v| <= 19/6`
guaranteed by

```text
s,t <= (19/6)^2 * 30 * 131 / (4*248) = 236455/5952.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n131_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_131`,

```text
chi(y) - chi(x) >= (496/131) g_R,
g_R =
  s0 - t1
  - (248/(6*50*131)) s1^2
  + (1/432)(496^2/1800) s0^3/131^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  1 - (s1+t1)/131
  + (1/12)(496/50)(s0^2+t0^2)/131^2
  - (1/360)(496^2/1800)(s1^3+t1^3)/131^3.
```

Moreover every retained cell has `g_R > 0`, `h_R > 0`, and
`496^2 h_R^2 - 4 > 0`.

Proof.  Lemmas 2 and 3 give the displayed formulas.  The endpoint choices
use lower endpoints in positive terms and upper endpoints in negative
terms.  The exact rational checker verifies the positivity comparisons over
all retained cells.

## Proposition 6 (rectangular radial mass lower)

For a retained cell `R=[s0,s1] x [t0,t1]`, define

```text
F_R =
  g_R^2 * (11/4)^(-ceil(s1+t1))
  * (s1^124-s0^124)/124
  * (t1^124-t0^124)/124
  * 1/(123!)^2
  / 124.
```

This is a lower bound for the normalized radial mass factor on `R`.

Proof.  On the cell the gap factor is at least `g_R^2`, and
`exp(-(s+t)) >= (11/4)^(-ceil(s1+t1))` by `e < 11/4`.  Since
`a=124`, `Gamma(124)=123!`, giving the displayed factor.

## Lemma 7 (exact-rational arithmetic check)

Let

```text
A0_lower = A0_prefactor(E8) * (113/355)^8,
Sine_R = (11/4)^(-ceil((4*248/(21*131))(s1+t1))),
Neg(131) = 520 * 16^131.
```

The rational lower ratio

```text
R131 =
  sum_R A0_lower * 496^131 * 131^(-250)
        * F_R * h_R^131
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(131)
```

satisfies `R131 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n131_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output
includes

```text
root_identity_constants: OK
retained_cells=17797
log10_ratio131_lower=1.59
log10_best_cell_ratio131=-1.41
positive_gaps: OK
positive_h131: OK
positive_direct131: OK
ratio131_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=131))

The rectangular-union positive contribution to `D_E8(131)` dominates the
negative-region bound `520*16^131`.  Hence `D_E8(131) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=131`.  Lemma 7 checks that this lower bound is larger than the
negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=131)

Assume `E8-FiniteBridge(131)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 131`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(133) >= 4 Q_3^{E8}(131)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_TAIL.md` gives
`Q_3^{E8}(n+2) >= 4 Q_3^{E8}(n)` for every odd `n >= 133`.  Induction over
odd indices from the assumed finite bridge gives the claim.
