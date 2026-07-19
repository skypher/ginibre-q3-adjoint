# DCT-RectTrig(E8, n=91) certificate

This note records the rectangular certificate for

```text
D_E8(91) = Q_3^{E8}(93) - 4 Q_3^{E8}(91).
```

## Lemma 1 (E8 constants and scalar bounds)

Use the E8 constants

```text
d = 248, |R_+| = 120, r = 8, alpha = 250, a = d/2 = 124.
```

On `0 <= w <= (13/3)^2`, set

```text
L_91(w) =
  (99119225/100000000)w
  - (7802858/100000000)w^2
  + (1562/890661)w^3,

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
L_91(w) <= 2(1-cos sqrt(w)) <= U(w),
C(w) <= 2 cos sqrt(w).
```

Proof.  The checker `character_ring_iter/direct_chain_rect_e8_n91_certificate.py`
verifies the three residual polynomials by exact Sturm root counts and returns

```text
n91_gap_lower_sturm: OK
n91_gap_upper_sturm: OK
n91_cos_lower_sturm: OK
root_identity_constants: OK
```

## Lemma 2 (sextic Weyl-sine bound)

For `|z| <= 13/3`,

```text
log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2880 - z^6/115000.
```

Proof.  With `y=(z/2)^2`, the checker verifies exact Sturm positivity for the
Taylor-side residual of

```text
sin(sqrt(y))/sqrt(y) * exp(y/6 + y^2/180 + 16y^3/28750) - 1
```

on `0 <= y <= (13/6)^2`.  It returns

```text
sine_sextic_13_3_sturm: OK
```

## Lemma 3 (optimized degree-27 negative majorant)

Let `S=X+Y`, where `X,Y` are independent E8 adjoint characters, and let

```text
J_k = int int (X-Y)^2 S^k dX dY.
```

There is a degree-27 polynomial

```text
P_27(S) = sum_{j=0}^{27} b_j T_j((2S-503)/489),
```

with the rational coefficients recorded in
`direct_chain_rect_e8_n91_certificate.py`, such that on the odd negative-sign
region

```text
S in [-16,-2] union [0,2],
```

one has

```text
|S^91(S^2-4)| <= 16^79 * S^12(S^2+4) * P_27(S)^2.
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
`references/oeis_A179663.txt`.  The exact integral is

```text
int int (X-Y)^2 S^12(S^2+4)P_27(S)^2 dX dY
= 154899922705202397910663859891726275260116389741415466199911998543539776953588818472524235571552289212256879365191595587107363092116739932995012713161412017164982537410758765023648674235489
  / 3865401811617780810470268062535894321947867605658317875555118084768343566134415279121319626970868866380880938847280490359795809438691655382933018750000000000000000000000000000000000000000.
```

Thus the negative part is bounded by `16^79` times this integral.  The checker
returns

```text
moment_log_0_70_oeis_overlap_0_30: OK
optimized_negative_moment_integral: OK
```

## Lemma 4 (rectangular union)

Let `U_91` be the union of `1/20 x 1/20` cells

```text
R(s0,t0) = [s0,s0+1/20] x [t0,t0+1/20],
s0 >= 40, t0 >= 30,
s0+1/20, t0+1/20 <= (13/3)^2 * 30 * 91 / (4*248).
```

Retain exactly those cells for which the gap, average-character lower bound,
and direct factor in Lemma 5 are positive.  The exact enumeration has `69502`
cells and cap

```text
(13/3)^2 * 30 * 91 / (4*248) = 76895/1488.
```

Proof.  This is the exact enumeration in
`direct_chain_rect_e8_n91_certificate.py`, which returns

```text
cell_count_69502: OK
all_cells_in_cap: OK
```

## Lemma 5 (local lower bounds on retained cells)

For each retained cell `R=[s0,s1] x [t0,t1]`, define `g_R` by the lower
polynomial `L_91` at `s0` minus the upper polynomial `U` at `t1`, using the
E8 quartic and sextic root identities.  Define `h_R` by the polynomial `C`,
using lower endpoints in positive terms and upper endpoints in negative
terms.  Then every retained cell satisfies

```text
g_R > 0, h_R > 0, 496^2 h_R^2 - 4 > 0.
```

Proof.  The checker verifies these exact rational inequalities over all
retained cells and returns

```text
positive_gaps: OK
positive_h91: OK
positive_direct91: OK
```

The minima are

```text
gap91_min = 29512424195754583/2676615873750000000000,
h91_min = 1173568851885107429/4171553750000000000,
direct_factor91_min =
  1376980911433068867147103318870990041
  / 70734670465087890625000000000000.
```

## Proposition 6 (positive rectangular lower bound)

For every retained cell, combine Lemmas 1, 2, and 5 with the fractional-part
exponential lower bound using `e < 1457/536`.  Summing the resulting exact
rational lower bounds over `U_91` gives

```text
log10_positive91_lower = 96.86.
```

The union radial fraction lower bound and largest sine exponent are

```text
log10_union_fraction_lower = -38.10,
sine_exponent_max = 16853426542806173/348190171875000.
```

Proof.  The checker performs the exact rational summation and returns

```text
e_upper_1457_536: OK
positive_radial_exp_lower: OK
positive_sine_exp_lower: OK
```

## Theorem 7 (E8 bridge step at n=91)

The direct Chain bridge step satisfies

```text
D_E8(91) = Q_3^{E8}(93) - 4 Q_3^{E8}(91) >= 0.
```

Proof.  Lemma 3 gives

```text
log10_negative91_upper = 96.73.
```

Proposition 6 gives

```text
log10_positive91_lower = 96.86.
```

The exact rational ratio is greater than `1`; the checker returns

```text
log10_ratio91_lower = 0.13
ratio91_gt_1: OK
```
