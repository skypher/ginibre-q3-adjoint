# DCT-RectTrig(E8, n=89) certificate

This note records the rectangular certificate for

```text
D_E8(89) = Q_3^{E8}(91) - 4 Q_3^{E8}(89).
```

## Lemma 1 (E8 constants and scalar bounds)

Use

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (14/3)^2`, set

```text
L_89(w) =
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
L_89(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n89_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n89_gap_lower_sturm: OK
n89_gap_upper_sturm: OK
n89_cos_lower_sturm: OK
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

## Lemma 3 (degree-28 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-28 polynomial

```text
P_28(S) = sum_{j=0}^{28} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n89_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^89(S^2-4)| <= 16^77 * S^12(S^2+4) * P_28(S)^2.
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
`references/oeis_A179663.txt`.  The two new values `m_71,m_72` are recorded
in `references/arxiv_2412_21189_e8_m71_m72.txt` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189.  The exact integral is

```text
int int (X-Y)^2 S^12(S^2+4)P_28(S)^2 dX dY
= 110574439399728242778064347972737777648628809018978860670416750279550250166056988167096970632724370157057840918236970893796884473467113003625506475576422917439930098480783344439652958219223
  / 62390165395220237149681115433261077065698483356952027438641538929482648026739857327216912260160826558355727591021787674134419528043283644560507030415015625000000000000000000000000000000000.
```

Thus

```text
|D_-(89)| <= 16^77 * int int (X-Y)^2 S^12(S^2+4)P_28(S)^2 dX dY.
```

## Lemma 4 (positive rectangular region)

Let

```text
delta = 14/3, grid_step = 1/20,
S_cap = delta^2 * 30*89/(4*248) = 21805/372.
```

Use rectangles with

```text
s0 >= 40, t0 >= 30, s1=s0+1/20, t1=t0+1/20,
s1,t1 <= S_cap,
```

retaining exactly those cells whose local gap, local average-character lower
bound, and direct factor are positive.  The checker retains `130740` cells.
For each retained cell the positive contribution is bounded below using
Lemmas 1 and 2, the exact E8 root identities, `e < 1457/536`, and the
Macdonald-Mehta `A_0` lower bound.

Proof.  The checker returns

```text
retained_cells=130740 cap=21805/372
log10_union_fraction_lower=-31.17
gap89_min=4441364733419/44555625000000000
h89_min=9192686056706299193/40976323125000000000
direct_factor89_min=84478176939489990409318539166032451249/6824999418104553222656250000000000
sine_exponent_max=14329109831533238/257258609296875
positive_gaps: OK
positive_h89: OK
positive_direct89: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Proposition 5 (DCT-RectTrig(E8, n=89))

The rectangular positive contribution dominates the negative-region bound:

```text
log10_positive89_lower=93.44
log10_negative89_upper=92.97
log10_ratio89_lower=0.47
ratio89_gt_1: OK
```

Proof.  Lemma 4 supplies the positive lower bound, and Lemma 3 supplies the
negative upper bound.  The exact rational comparison in the checker is
`positive_total / negative_bound > 1`.

Therefore

```text
D_E8(89) >= 0.
```

## Corollary 6 (conditional E8 odd tail from n=89)

Assume `E8-FiniteBridge(89)`: `Q_3^{E8}(n) >= 0` for every odd
`67 < n <= 89`.  Then `Q_3^{E8}(n) >= 0` for every odd `n >= 89`.

Proof.  Proposition 5 gives `D_E8(89) >= 0`, and the previously recorded
bridge certificates give the direct Chain steps for every odd
`n=91,93,...,131`, followed by the E8 rectangular tail from odd `n >= 133`.
Thus the direct Chain propagates nonnegativity from the finite bridge window
to every odd `n >= 89`.
