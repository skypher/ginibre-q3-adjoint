# `SU(2)_13` and the rank-seven odd-orbit frontier

Date: 2026-07-25

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

This census uses long-double spectral arithmetic and is a frontier diagnostic.
Every theorem below is instead replayed with rational Sturm isolation, exact
rational capacities, directed-rounding MPFR logarithms where needed, and
integral fusion arithmetic for finite leaves.

## Exact frontier: 175 closed failure keys

The exact program now closes 175 direct-Hall failure keys.

### Initial 172-key package

The earlier exact package proves:

```text
first ray                         1
second ray                        1
complete support-164 chamber      3
support-167 ray                   1
support-169 regimes               3
generic multivariable replay    163
                                ---
total                           172.
```

The first two rays are

```text
(B_1^-)^(1+2p) B_4^-,
(B_1^-)^(2+2p) B_5^+,
```

and both are proved by exact finite leaves followed by denominator-100
weighted AM-GM tails. The package also proves the full two-variable
`B_1^-B_5^-` chamber, decorated rays, and 163 multivariable regimes through
supports `178,187,196,205,214,223,232,241,250`.

Replay the initial package with

```text
character_ring_iter/replay_su2_o13_initial_frontier_package.sh
```

Its SHA-256 is

```text
a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2.
```

### Support 245: complete

Support code `245` has signs

```text
B_1 minus, B_6 plus,
```

with minimum exponents two and one. All three non-Hall residual faces are now
proved:

```text
(B_1^-)^(4+2p) B_6^+,
(B_1^-)^2 (B_6^+)^(3+2q),
(B_1^-)^(4+2p) (B_6^+)^(3+2q),          p,q>=0.
```

Exact transfer arithmetic found:

```text
B_1 ray zeros: p=0 only;
B_6 ray zeros: q=0,1,2 only;
121 x 121 plane scan: only (0,0) is zero.
```

The all-exponent proof uses a slope-three decomposition of the residual plane:

```text
2 cones,
12 diagonal strips,
2 column tails,
14 finite column points,
```

plus the two standalone ray tails. In total the strict verifier checks

```text
18 denominator-100 AM-GM pieces,
20 exact finite leaves.
```

It reports

```text
SU2_O13_SUPPORT245_EXACT PASS
pieces=18 denominator=100 leaves=20 zeros=5
minimum_leaf=0 slope=3 strips=12.
```

The detailed theorem is

```text
SU2_O13_SUPPORT245_GKS2_2026_07_25.md.
```

Replay with

```text
character_ring_iter/replay_su2_o13_support245_package.sh
```

The deterministic package SHA-256 is

```text
da31783b2297c692d7dbc4b10e9f41bc7a03b07b664b43c137403ad06995c76f.
```

## Next precise target

The next support in exact support order is `248`. The diagnostic failure keys
are

```text
support=248 parity=4 residual=5
support=248 parity=4 residual=7
support=248 parity=6 residual=1
support=248 parity=6 residual=4
support=248 parity=6 residual=5
support=248 parity=6 residual=7.
```

Their residual variables involve `B_1`, `B_2`, and `B_6`. Full `O_13`, and
hence full `SU(2)_13`, remains open.
