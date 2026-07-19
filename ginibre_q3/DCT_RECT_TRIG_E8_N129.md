# DCT-RectTrig(E8, n=129) certificate

This note records the fixed wide-cap rectangular certificate for the single
E8 bridge step

```text
D_E8(129) = Q_3^{E8}(131) - 4 Q_3^{E8}(129).
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
by `direct_chain_rect_e8_n129_certificate.py`.

## Lemma 2 (wide-cap trigonometric bounds)

For `|z| <= 19/6`,

```text
2(1-cos z) >= z^2 - z^4/12 + z^6/432,
2 cos z >= 2 - z^2 + z^4/12 - z^6/360,
sin(z/2)/(z/2) >= exp(-z^2/21).
```

Proof.  This is Lemma 2 of `DCT_RECT_TRIG_E8_N131.md`.

## Lemma 3 (rational exponential lower bound)

For `x >= 0`, write `x=k+r`, where `k=floor(x)` and `0 <= r < 1`.  Then

```text
exp(-x) >= (1457/536)^(-k) * (1 - r + r^2/2 - r^3/6).
```

Proof.  The checker verifies

```text
sum_{j=0}^7 1/j! + 1/(7!*7) < 1457/536.
```

Since the Taylor tail after degree `7` is bounded by `1/(7!*7)`, this gives
`e < 1457/536`.  On `0 <= r < 1`, the alternating Taylor expansion gives

```text
exp(-r) >= 1 - r + r^2/2 - r^3/6.
```

Multiplying the two lower bounds gives the displayed inequality.

## Lemma 4 (wide-cap rectangular union)

Let `U_129` be the union of `1/40 x 1/40` cells

```text
R(s0,t0) = [s0,s0+1/40] x [t0,t0+1/40],
s0 in {37, 1481/40, ...},
t0 in {26, 1041/40, ..., 1239/40},
```

retained by the positive rational gap in Lemma 5.  The exact enumeration has
`12719` cells.  Each cell lies inside the root cap `|alpha.v| <= 19/6`
guaranteed by

```text
s,t <= (19/6)^2 * 30 * 129 / (4*248) = 77615/1984.
```

Proof.  This is the exact enumeration performed by
`direct_chain_rect_e8_n129_certificate.py`.

## Lemma 5 (local character lower bounds on the union)

For a retained cell `R=[s0,s1] x [t0,t1]` in `U_129`,

```text
chi(y) - chi(x) >= (496/129) g_R,
g_R =
  s0 - t1
  - (248/(6*50*129)) s1^2
  + (1/432)(496^2/1800) s0^3/129^2,
```

and

```text
(chi(x)+chi(y))/496 >= h_R,
h_R =
  1 - (s1+t1)/129
  + (1/12)(496/50)(s0^2+t0^2)/129^2
  - (1/360)(496^2/1800)(s1^3+t1^3)/129^3.
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
Sine_R = ExpLower((4*248/(21*129))(s1+t1)),
Neg(129) = 520 * 16^129.
```

The rational lower ratio

```text
R129 =
  sum_R A0_lower * 496^129 * 129^(-250)
        * F_R * h_R^129
        * (496^2 h_R^2 - 4)
        * Sine_R
  / Neg(129)
```

satisfies `R129 > 1`.

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_n129_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output
includes

```text
root_identity_constants: OK
e_upper_1457_536: OK
retained_cells=12719
log10_ratio129_lower=0.35
log10_best_cell_ratio129=-2.90
positive_gaps: OK
positive_h129: OK
positive_direct129: OK
ratio129_gt_1: OK
```

## Theorem 8 (DCT-RectTrig(E8, n=129))

The rectangular-union positive contribution to `D_E8(129)` dominates the
negative-region bound `520*16^129`.  Hence `D_E8(129) >= 0`.

Proof.  Lemmas 4-6 give the positive rectangular-union lower bound at
`n=129`.  Lemma 7 checks that this lower bound is larger than the
negative-region bound.

## Corollary 9 (conditional E8 odd tail from n=129)

Assume `E8-FiniteBridge(129)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 129`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives `Q_3^{E8}(131) >= 4 Q_3^{E8}(129)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_N131.md` gives
`Q_3^{E8}(133) >= 4 Q_3^{E8}(131)`.  Theorem 8 of
`DCT_RECT_TRIG_E8_TAIL.md` gives
`Q_3^{E8}(n+2) >= 4 Q_3^{E8}(n)` for every odd `n >= 133`.  Induction over
odd indices from the assumed finite bridge gives the claim.
