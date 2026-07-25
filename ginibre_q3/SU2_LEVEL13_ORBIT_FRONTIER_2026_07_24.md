# `SU(2)_13` and the rank-seven odd-orbit frontier

Date: 2026-07-24

## Exact reduction

By odd-level simple-current lifting,

```text
full SU(2)_13 scalar GKS2*  <=>  scalar GKS2* in the rank-seven orbit ring O_13.
```

Use the even-lift orbit basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6,
B_4=V_8, B_5=V_10, B_6=V_12.
```

## Parallel direct-transport census

A five-thread C++ diagnostic enumerated the support/sign/parity chambers and
all residual coordinate faces. After excluding the already proved
separated-pole sector, it found

```text
admissible chambers       6,783
residual regimes        243,008
direct Hall regimes     208,084
non-Hall regimes         34,924
separated-pole regimes        4.
```

This census uses long-double spectral arithmetic and is a frontier diagnostic,
not the exact theorem below.

## First direct-Hall obstruction: proved

The first non-Hall residual ray has

```text
support=56, parity=3, residual=1,
```

which decodes to

```text
F_p=[B_0 tensor B_0]
    (B_1 tensor 1-1 tensor B_1)^(1+2p)
    (B_4 tensor 1-1 tensor B_4),          p>=0.
```

**Theorem 1.** `F_p>=0` for every `p>=0`.

The exact verifier uses the seven real roots of

```text
x^7-6x^6+9x^5+5x^4-15x^3+5x
 =x(x^2-x-1)(x^4-5x^3+5x^2+5x-5).
```

At a trace node,

```text
B_1=x,
B_4=x^4-3x^3+3x,
weight=(3-x)/15.
```

Rational Sturm sequences isolate all seven roots. One nominal spectral pair
vanishes identically: the two roots of `x^2-x-1` both have `B_4=-1`.
The remaining spectral sum has ten positive and ten negative terms.

Exact integral orbit-fusion arithmetic gives

```text
F_0=0,
F_1=0,
F_2=8.
```

After shifting the residual floor to `p=2`, a denominator-100 capacitated
weighted AM-GM allocation proves the entire tail. Every coefficient,
squared base, geometric inequality, and capacity inequality is checked with
rational interval arithmetic. Thus the two zero leaves and the AM-GM tail
cover every `p>=0`.

The strict replay reports

```text
SU2_O13_FIRST_RAY_EXACT PASS
roots=7 spectral_pairs=21 positives=10 negatives=10
floor=2 denominator=100 leaves=2 leaf_values=0,0 tail_first=8.
```

A separate exact integer scan through `p=2000` found no additional zero and no
negative value. The relevant scalar sequence has a degree-20 minimal linear
recurrence; its first values are

```text
0, 0, 8, 196, 3456, 54430, 816270, 11950512, ...
```

The recurrence is a regression and structural diagnostic; the AM-GM
certificate is the all-exponent proof.

## Next precise target

After removing the separated-pole sector and the proved ray above, the next
direct-Hall obstruction is

```text
support=83, parity=2, residual=1.
```

It decodes to the one-variable tail

```text
(B_1^-)^(2+2r) B_5^+,       r>=1,
```

with the `r=0` face handled separately. This is the next exact rank-seven
orbit target. Full `O_13`, and hence full `SU(2)_13`, remains open.
