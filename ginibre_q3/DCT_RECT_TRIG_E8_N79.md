# DCT-RectTrig(E8, n=79) certificate

This note records the rectangular certificate for

```text
D_E8(79) = Q_3^{E8}(81) - 4 Q_3^{E8}(79).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, use the same cap-adjusted scalar polynomials as in
`DCT_RECT_TRIG_E8_N89.md`:

```text
L_79(w) =
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
L_79(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n79_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n79_gap_lower_sturm: OK
n79_gap_upper_sturm: OK
n79_cos_lower_sturm: OK
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

## Lemma 3 (degree-37 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-37 polynomial

```text
P_37(S) = sum_{j=0}^{37} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n79_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^79(S^2-4)| <= 16^57 * S^22(S^2+4) * P_37(S)^2.
```

Proof.  After division by the positive factor `x^22` on `S=-x in [-16,-2]`,
the negative interval residual is checked by exact Bernstein coefficient
positivity on `[2,16]`.  After division by `S^22` on `S in [0,2]`, the
positive interval residual is checked by exact Bernstein coefficient
positivity on `[0,2]`.  The checker returns

```text
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
```

The moment source check uses `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log`, with `m_0,...,m_30` checked against
`references/oeis_A179663.txt`.  The thirty new values `m_71,...,m_100` are
recorded in `references/arxiv_2412_21189_e8_m71_m100.txt` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189.  The exact integral is

```text
int int (X-Y)^2 S^22(S^2+4)P_37(S)^2 dX dY
= 29822698848482859512999327002571175457028680982472229890974727314690982554293289858766526952918578162868771228943226337357855616889719317471018088299258706043797188858166167899916127307423211905059325604991005365894898566191400370099111060755564260674533
  / 10206000262066657917098699194004437811274933703122998614336992087278816674295977124331635069516810245691721021685123596041961334507464344450115000019979980691494330192138861070759341811883181826357841000000000000000000000000000000000000000000000000.
```

Thus

```text
|D_-(79)| <= 16^57 * int int (X-Y)^2 S^22(S^2+4)P_37(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*79/(4*248) = 19355/372.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `68095` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=68095 cap=19355/372
log10_union_fraction_lower=-38.71
gap79_min=2975877287857/35105625000000000
h79_min=2141716707210133747/9552630625000000000
direct_factor79_min=4585466765860003846337353449628260009/370922020753479003906250000000000
sine_exponent_max=3340401796495402/59973572109375
positive_gaps: OK
positive_h79: OK
positive_direct79: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=79))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive79_lower=75.89
log10_negative79_upper=75.10
log10_ratio79_lower=0.79
ratio79_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(79) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=79)

Assume `E8-FiniteBridge(79)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 79`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 79`.

Proof.  Proposition 5 gives `D_E8(79) >= 0`; hence `Q_3^{E8}(81) >= 0`
follows from `Q_3^{E8}(79) >= 0`.  Corollary 6 of
`DCT_RECT_TRIG_E8_N81.md` then propagates nonnegativity from the finite bridge
window to every odd `n >= 81`.
