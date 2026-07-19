# DCT-RectTrig(E8, n=77) certificate

This note records the rectangular certificate for

```text
D_E8(77) = Q_3^{E8}(79) - 4 Q_3^{E8}(77).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, use the same cap-adjusted scalar polynomials as in
`DCT_RECT_TRIG_E8_N89.md`:

```text
L_77(w) =
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
L_77(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n77_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n77_gap_lower_sturm: OK
n77_gap_upper_sturm: OK
n77_cos_lower_sturm: OK
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
`direct_chain_rect_e8_n77_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^77(S^2-4)| <= 16^55 * S^22(S^2+4) * P_37(S)^2.
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
= 275850426105466974599271222814594639961288531779232471683078843699472794505799470453340725385828976397284323226218809282233398263426608337113956109325735037174032212343506831322288989678916089949771483864475633877431448353313163971042082653208715595538053
  / 255150006551666447927467479850110945281873342578074965358424802181970416857399428108290876737920256142293025542128089901049033362686608611252875000499499517287358254803471526768983545297079545658946025000000000000000000000000000000000000000000000000.
```

Thus

```text
|D_-(77)| <= 16^55 * int int (X-Y)^2 S^22(S^2+4)P_37(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*77/(4*248) = 18865/372.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `57817` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=57817 cap=18865/372
log10_union_fraction_lower=-40.40
gap77_min=1815834031/275625000000000
h77_min=3964993922572267201/17690653750000000000
direct_factor77_min=15716088369021122430214894393340374401/1272109253472900390625000000000000
sine_exponent_max=2749674766803241/49362630625000
positive_gaps: OK
positive_h77: OK
positive_direct77: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=77))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive77_lower=72.46
log10_negative77_upper=72.26
log10_ratio77_lower=0.19
ratio77_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(77) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=77)

Assume `E8-FiniteBridge(77)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 77`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 77`.

Proof.  Proposition 5 gives `D_E8(77) >= 0`; hence `Q_3^{E8}(79) >= 0`
follows from `Q_3^{E8}(77) >= 0`.  Corollary 6 of
`DCT_RECT_TRIG_E8_N79.md` then propagates nonnegativity from the finite bridge
window to every odd `n >= 79`.
