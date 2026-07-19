# DCT-RectTrig(E8, n=81) certificate

This note records the rectangular certificate for

```text
D_E8(81) = Q_3^{E8}(83) - 4 Q_3^{E8}(81).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, use the same cap-adjusted scalar polynomials as in
`DCT_RECT_TRIG_E8_N89.md`:

```text
L_81(w) =
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
L_81(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n81_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n81_gap_lower_sturm: OK
n81_gap_upper_sturm: OK
n81_cos_lower_sturm: OK
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

## Lemma 3 (degree-42 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-42 polynomial

```text
P_42(S) = sum_{j=0}^{42} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n81_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^81(S^2-4)| <= 16^69 * S^12(S^2+4) * P_42(S)^2.
```

Proof.  After division by the positive factor `x^12` on `S=-x in [-16,-2]`,
the negative interval residual is checked by exact Bernstein coefficient
positivity on the 8 equal subintervals of `[2,16]`.  After division by `S^12`
on `S in [0,2]`, the positive interval residual is checked by exact Bernstein
coefficient positivity on `[0,2]`.  The checker returns

```text
optimized_negative_interval_bernstein_8: OK
optimized_positive_interval_bernstein: OK
```

The moment source check uses `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log`, with `m_0,...,m_30` checked against
`references/oeis_A179663.txt`.  The thirty new values `m_71,...,m_100` are
recorded in `references/arxiv_2412_21189_e8_m71_m100.txt` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189.  The exact integral is

```text
int int (X-Y)^2 S^12(S^2+4)P_42(S)^2 dX dY
= 162988295159444920217175613977736273000922799417640375650345408112354141693199993604285605069801764134398766069421665451134117097220188593300942224975121288698309839236172476378298497541672298842570954921593146535995554007868966392982351229089785126310293950485983753
  / 1662275213223649722955032911709714874544854221093257063513796535114453366206250958407904534456796140022503212421180405346282041850376672598241125169818513010264936158883955387393292011765444675330836543152256442397205420466966875000000000000000000000000000000000000000000.
```

Thus

```text
|D_-(81)| <= 16^69 * int int (X-Y)^2 S^12(S^2+4)P_42(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*81/(4*248) = 6615/124.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `79060` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=79060 cap=6615/124
log10_union_fraction_lower=-37.10
gap81_min=2237512207/164025000000000
h81_min=1733570508555555541/7722502031250000000
direct_factor81_min=3004297063271076240979522552635802681/242411215622806549072265625000000
sine_exponent_max=10795237902751954/193934446171875
positive_gaps: OK
positive_h81: OK
positive_direct81: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=81))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive81_lower=79.34
log10_negative81_upper=79.08
log10_ratio81_lower=0.27
ratio81_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(81) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=81)

Assume `E8-FiniteBridge(81)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 81`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 81`.

Proof.  Proposition 5 gives `D_E8(81) >= 0`, and the previously recorded
bridge certificates give the direct Chain steps for every odd
`n=83,85,...,131`, followed by the E8 rectangular tail from odd `n >= 133`.
Thus the direct Chain propagates nonnegativity from the finite bridge window
to every odd `n >= 81`.
