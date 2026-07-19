# DCT-RectTrig(E8, n=85) certificate

This note records the rectangular certificate for

```text
D_E8(85) = Q_3^{E8}(87) - 4 Q_3^{E8}(85).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, use the same cap-adjusted scalar polynomials as in
`DCT_RECT_TRIG_E8_N89.md`:

```text
L_85(w) =
  (99119225/100000000)w
  - (7802858/100000000)w^2
  + (169078/100000000)w^3,

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
L_85(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n85_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n85_gap_lower_sturm: OK
n85_gap_upper_sturm: OK
n85_cos_lower_sturm: OK
root_identity_constants: OK
```

## Lemma 2 (sextic Weyl-sine bound)

For `|z| <= 14/3`,

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2880 - z^6/103800.
```

Proof.  With `y=(z/2)^2`, the checker verifies exact Sturm positivity for the
Taylor-side residual of

```text
sin(sqrt(y))/sqrt(y) * exp(y/6 + y^2/180 + 8y^3/12975) - 1
```

on `0 <= y <= (7/3)^2`.  It returns

```text
sine_sextic_14_3_sturm: OK
```

## Lemma 3 (degree-34 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-34 polynomial

```text
P_34(S) = sum_{j=0}^{34} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n85_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^85(S^2-4)| <= 16^73 * S^12(S^2+4) * P_34(S)^2.
```

Proof.  After division by the positive factor `x^12` on `S=-x in [-16,-2]`,
the negative interval residual is checked by exact Bernstein coefficient
positivity on `[2,16]`.  After division by `S^12` on `S in [0,2]`, the
positive interval residual is checked by exact Bernstein coefficient
positivity on each of the 32 equal subintervals of `[0,2]`.  The checker
returns

```text
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein_32: OK
```

The moment source check uses `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log`, with `m_0,...,m_30` checked against
`references/oeis_A179663.txt`.  The fourteen new values `m_71,...,m_84` are
recorded in `references/arxiv_2412_21189_e8_m71_m84.txt` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189.  The exact integral is

```text
int int (X-Y)^2 S^12(S^2+4)P_34(S)^2 dX dY
= 65148687217721097618389819413073654708609879754137122790264611448139731858033630148113111544815479392877378205889026723166072799733022265663747429139658760314718523759506158358410320170203871509898553101174638646758322907219
  / 4665330740770604793149491012159339221112632191672044571621781307705608579427809365087496794344569859695546279931419517931942002264051301274362377419277716193612079960316166376689253006250000000000000000000000000000000000000000.
```

Thus

```text
|D_-(85)| <= 16^73 * int int (X-Y)^2 S^12(S^2+4)P_34(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*85/(4*248) = 20825/372.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `103495` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=103495 cap=20825/372
log10_union_fraction_lower=-34.01
gap85_min=3855342619/37353515625000
h85_min=213613375085470649/951893750000000000
direct_factor85_min=45615941612451871037465981840481201/3683100738525390625000000000000
sine_exponent_max=147898097562783/2656090625000
positive_gaps: OK
positive_h85: OK
positive_direct85: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=85))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive85_lower=86.34
log10_negative85_upper=86.05
log10_ratio85_lower=0.29
ratio85_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(85) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=85)

Assume `E8-FiniteBridge(85)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 85`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 85`.

Proof.  Proposition 5 gives `D_E8(85) >= 0`, and the previously recorded
bridge certificates give the direct Chain steps for every odd
`n=87,89,...,131`, followed by the E8 rectangular tail from odd `n >= 133`.
Thus the direct Chain propagates nonnegativity from the finite bridge window
to every odd `n >= 85`.
