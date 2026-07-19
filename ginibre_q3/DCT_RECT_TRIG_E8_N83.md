# DCT-RectTrig(E8, n=83) certificate

This note records the rectangular certificate for

```text
D_E8(83) = Q_3^{E8}(85) - 4 Q_3^{E8}(83).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, use the same cap-adjusted scalar polynomials as in
`DCT_RECT_TRIG_E8_N89.md`:

```text
L_83(w) =
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
L_83(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n83_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n83_gap_lower_sturm: OK
n83_gap_upper_sturm: OK
n83_cos_lower_sturm: OK
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

## Lemma 3 (degree-38 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-38 polynomial

```text
P_38(S) = sum_{j=0}^{38} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n83_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^83(S^2-4)| <= 16^71 * S^12(S^2+4) * P_38(S)^2.
```

Proof.  After division by the positive factor `x^12` on `S=-x in [-16,-2]`
and by `S^12` on `S in [0,2]`, the two residual polynomials are checked by
exact Bernstein coefficient positivity.  The checker returns

```text
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
```

The moment source check uses `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log`, with `m_0,...,m_30` checked against
`references/oeis_A179663.txt`.  The twenty-two new values `m_71,...,m_92`
are recorded in `references/arxiv_2412_21189_e8_m71_m92.txt` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189.  The exact integral is

```text
int int (X-Y)^2 S^12(S^2+4)P_38(S)^2 dX dY
= 69506923908375742396638750904720698654253961502478284262296998242926162620300793054819096177856546091509678268326171062873371758355860924644968831631635906940012075850490819786026964402425828891779371420933273301225535316551879001356040185988627
  / 61011724716641032694863951249238379346746835550611863791471897122554948049358208648682622736448229569001250560659210985228745906718984537731398722994440824073270393246860914952526114334982958037512832444025000000000000000000000000000000000000000000.
```

Thus

```text
|D_-(83)| <= 16^71 * int int (X-Y)^2 S^12(S^2+4)P_38(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*83/(4*248) = 20335/372.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `91138` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=91138 cap=20335/372
log10_union_fraction_lower=-35.50
gap83_min=1339382274583/19375312500000000
h83_min=119188825916243597/531761910000000000
direct_factor83_min=14201378624292119757409379043498409/1149399750125390625000000000000
sine_exponent_max=3720113027348407/66770426925000
positive_gaps: OK
positive_h83: OK
positive_direct83: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=83))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive83_lower=82.83
log10_negative83_upper=82.55
log10_ratio83_lower=0.28
ratio83_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(83) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=83)

Assume `E8-FiniteBridge(83)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 83`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 83`.

Proof.  Proposition 5 gives `D_E8(83) >= 0`, and the previously recorded
bridge certificates give the direct Chain steps for every odd
`n=85,87,...,131`, followed by the E8 rectangular tail from odd `n >= 133`.
Thus the direct Chain propagates nonnegativity from the finite bridge window
to every odd `n >= 83`.
