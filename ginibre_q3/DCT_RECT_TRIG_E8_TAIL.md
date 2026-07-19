# DCT-RectTrig(E8-tail) certificate

This note records the fixed rectangular certificate for the E8 direct Chain
tail above the current exact prefix.  Write

```text
D_E8(n) = Q_3^{E8}(n+2) - 4 Q_3^{E8}(n).
```

The certificate starts at odd `n0 = 133`.  It closes the analytic tail
subproblem from that point onward; the finite E8 bridge below `133` remains
a separate named obligation.

## Lemma 1 (E8 constants)

For E8 in the adjoint normalization used by the direct Chain ledger,

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, C = 8, a = d/2 = 124.
```

The direct-tail union below is tested at `n0 = 133`.

Proof.  These are the E8 entries in the direct-tail root-data table.  They
are replayed by `direct_chain_rect_e8_tail_certificate.py` from the same
`ROOT_DATA` table used by the A0 and truncated-mass checkers.

## Lemma 2 (elementary trigonometric bounds)

For `|z| <= pi`,

```text
2(1-cos z) <= z^2,
2(1-cos z) >= z^2 - z^4/12,
2 cos z >= 2 - z^2 + z^4/18,
sin(z/2)/(z/2) >= exp(-z^2/18).
```

Also,

```text
sin(z/2)/(z/2) >= exp(-z^2/21).
```

Proof.  The first four inequalities are Lemma 2 of
`DCT_RECT_TRIG_E7.md`; those statements are group-independent.  For the
sharper sine bound, put `x=z/2` and

```text
F(x) = log(sin x/x) + 4x^2/21.
```

Using

```text
cot x = 1/x - 2x sum_{m>=1} 1/(m^2 pi^2 - x^2),
```

we have

```text
F'(x)/x = 8/21 - 2 sum_{m>=1} 1/(m^2 pi^2 - x^2).
```

The sum is increasing on `0 <= x <= pi/2`, so `F'` changes sign at most
once.  Hence the minimum of `F` occurs at an endpoint.  At `0`, `F=0`.
At `pi/2`, the needed inequality is

```text
log(pi/2) < pi^2/21.
```

Using `pi < 22/7`, `x=(pi-2)/2 < 4/7`, and the alternating logarithm
series,

```text
log(pi/2) = log(1+x)
 <= 4/7 - (4/7)^2/2 + (4/7)^3/3 - (4/7)^4/4 + (4/7)^5/5
 < (333/106)^2/21 < pi^2/21.
```

The rational bound `e < 11/4` used below follows from

```text
e = sum_{k>=0} 1/k!
  < 1 + 1 + 1/2 + 1/6 + (1/24) sum_{j>=0} 5^{-j}
  = 87/32 < 11/4.
```

## Lemma 3 (E8 quartic root identity)

For E8, with `Q(u)=sum_{alpha>0}(alpha.u)^2`,

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50.
```

Proof.  The left side is a Weyl-invariant quartic polynomial, and E8 has no
basic invariant of degree `4`; hence it is a scalar multiple of `Q^2`.  To
fix the scalar, evaluate at a root `beta`.  The exact E8 root enumeration in
`direct_chain_rect_e8_tail_certificate.py` gives

```text
sum_{alpha>0}(alpha.beta)^2 = 60,
sum_{alpha>0}(alpha.beta)^4 = 72,
```

so the scalar is `72/60^2 = 1/50`.

## Lemma 4 (rectangular union containment)

Let

```text
Delta = 1/40,
s0 in {37, 1481/40, ..., 793/20},
t0 in {26, 1041/40, ..., 1239/40},
R(s0,t0) = [s0,s0+Delta] x [t0,t0+Delta].
```

Let `U_133` be the union of the cells `R(s0,t0)` whose rational lower gap
in Lemma 5 is positive at `n=133`.  The exact enumeration has `13558`
cells.  For every `n >= 133`, every cell of `U_133` lies inside the root cap
`|alpha.v| <= pi` guaranteed by

```text
s,t <= pi^2 kappa_E8 n/(4d),     kappa_E8 = 30.
```

Proof.  At `n=133`, the cap is `1995 pi^2/496`.  Since
`pi > 333/106`, this is larger than

```text
221223555/5573056 > 1587/40.
```

The largest retained `s` endpoint is at most `1587/40`, and the largest
retained `t` endpoint is at most `31`.  The cap increases with `n`.

## Lemma 5 (local character lower bounds on each cell)

On a cell `R=[s0,s1] x [t0,t1]` of Lemma 4, with the larger radial variable
called `s` and the smaller one called `t`,

```text
chi(y) - chi(x) >= (2d/n) g,
g = s0 - t1 - (248/(300n))s1^2,
```

and

```text
(chi(x)+chi(y))/(2d) >= h(n),
h(n) = 1 - (s1+t1)/n
       + (1/18)(496/50)(s0^2+t0^2)/n^2.
```

For the retained cells at `n=133`,

```text
min h(133) = 176378749/318402000,
min(496^2 h(133)^2 - 4)
 = 29894609978051543461/396014975015625 > 0.
```

Proof.  Lemmas 2 and 3 give

```text
d - chi(x) >= Q(x)/n - Q(x)^2/(600 n^2),
d - chi(y) <= Q(y)/n.
```

With `Q=2ds`, this gives the displayed `g` after substituting the rectangle
endpoints.  The lower bound for `chi(x)+chi(y)` follows from
`2cos z >= 2-z^2+z^4/18` and Lemma 3.  Exact enumeration of the retained
cells gives the displayed minima.

## Proposition 6 (rectangular-union radial mass lower)

For a retained cell `R=[s0,s1] x [t0,t1]`, the normalized radial mass factor
on that cell is at least

```text
F_R =
  g_R^2 * (11/4)^(-ceil(s1+t1))
  * (s1^124 - s0^124)/124
  * (t1^124 - t0^124)/124
  * 1/(123!)^2
  / 124.
```

The union lower mass is the sum of these `F_R` over the `13558` retained
cells.

Proof.  On each cell the gap factor is at least `g_R^2`, and
`exp(-(s+t)) >= (11/4)^(-ceil(s1+t1))` by `e < 11/4`.  Since
`a=124`, `Gamma(a)=123!`, giving the displayed rational lower bound.  The
retained cells are disjoint, so their lower bounds add.

## Lemma 7 (exact-rational arithmetic check)

Let

```text
A0_lower = A0_prefactor(E8) * (113/355)^8,
Neg(n) = 520 * 16^n.
```

For each retained cell, put

```text
Sine_R = (11/4)^(-ceil(496(s1+t1)/(21n))).
```

The rational lower ratio

```text
R133 =
  sum_R A0_lower * 496^133 * 133^(-250)
        * F_R * h_R(133)^133
        * (496^2 h_R(133)^2 - 4)
        * Sine_R
  / Neg(133)
```

satisfies `R133 > 1`.  Moreover, with `h_min = 176378749/318402000`,

```text
((248/8)h_min)^2 * (133/135)^250 > 1.
```

Proof.  The checker
`character_ring_iter/direct_chain_rect_e8_tail_certificate.py` performs
these comparisons as exact integer comparisons.  Its current output
includes

```text
quartic_identity_constant: OK
retained_cells=13558
log10_ratio133_lower=0.02
log10_best_cell_ratio133=-3.25
log10_odd_step_ratio_lower=0.85
h_monotone_from_n133: OK
```

## Theorem 8 (DCT-RectTrig(E8-tail))

For every odd `n >= 133`, the rectangular-union positive contribution to
`D_E8(n)` dominates the negative-region bound `520*16^n`.

Proof.  Lemmas 4-6 give the stated positive rectangular-union lower bound
at each `n >= 133`.  Lemma 7 checks that every cell has
`h_R(n) >= h_R(133)` for `n >= 133`, and
`496^2 h_R(n)^2-4 >= 496^2 h_R(133)^2-4`.  Lemma 7 gives the base
comparison at `n=133`.  The same lemma gives an odd-step lower ratio larger
than `1`, so the lower comparison propagates from `n` to `n+2` for all odd
`n >= 133`.

## Corollary 9 (conditional E8 odd tail)

Assume `E8-FiniteBridge(133)`: `Q_3^{E8}(n) >= 0` for every odd
`67 <= n <= 133`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 67`.

Proof.  Theorem 8 gives
`D_E8(n)=Q_3^{E8}(n+2)-4Q_3^{E8}(n) >= 0` for every odd `n >= 133`.
Thus `Q_3^{E8}(n+2) >= 4Q_3^{E8}(n)`, and induction over odd indices from
the assumed finite bridge gives the claim.
