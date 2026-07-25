# `SU(2)_13` and the rank-seven odd-orbit frontier

Date: 2026-07-25

## Exact reduction

Odd-level simple-current lifting gives

```text
full SU(2)_13 scalar GKS2*  <=>  scalar GKS2* in O_13.
```

Use the even-lift orbit basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6,
B_4=V_8, B_5=V_10, B_6=V_12.
```

## Diagnostic census

A five-thread C++ census found

```text
admissible chambers       6,783
residual regimes        243,008
direct Hall regimes     208,084
non-Hall regimes         34,924
separated-pole regimes        4.
```

The census is a frontier diagnostic. The exact results below use rational
Sturm isolation, rational interval arithmetic, directed-rounding MPFR
logarithms where necessary, and integral fusion arithmetic for finite leaves.

## Exact frontier: 211 closed failure keys

The cumulative exact program now proves 211 direct-Hall failure keys:

```text
initial package through selected support 250       172
complete support 245                                 3
complete supports 248,249,251                       17
complete support 259                                19
                                                     ---
total                                               211.
```

### Initial 172-key package

This package proves the first two rays, the complete `B_1^-B_5^-` chamber,
decorated rays, and 163 multivariable regimes through supports
`178,187,196,205,214,223,232,241,250`.

Replay:

```text
character_ring_iter/replay_su2_o13_initial_frontier_package.sh
```

SHA-256:

```text
a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2.
```

### Support 245

All three support-245 faces are proved by two ray tails and a slope-three plane
decomposition. The strict replay checks 18 denominator-100 AM-GM pieces and
20 exact leaves.

See

```text
SU2_O13_SUPPORT245_GKS2_2026_07_25.md
character_ring_iter/replay_su2_o13_support245_package.sh.
```

### Supports 248, 249, and 251

All 17 regimes at these supports are proved. The mixed support-251 faces use
31 Farey-fan pieces and a four-piece translated-tail recursion.

See

```text
SU2_O13_SUPPORT248_251_GKS2_2026_07_25.md
character_ring_iter/replay_su2_o13_support248_251_package.sh.
```

Package SHA-256:

```text
0422d1fd71b20e3229730abfcd07f4b67d17497a1566e6f9c5b053733f37fdaf.
```

### Support 259

All 19 diagnostic failure regimes at support `259` have global rational
weighted-AM-GM allocations:

```text
18 denominator-100 tables,
1 denominator-200 table,
0 subdivisions,
0 finite leaves.
```

The strict exact verifier reports

```text
SU2_O13_SUPPORT259_EXACT PASS cases=19.
```

See

```text
SU2_O13_SUPPORT259_GKS2_2026_07_25.md
certificates/su2_o13_support259_exact.log.
```

The local source package was replayed but still needs repository restoration
because the connector altered larger text chunks during the attempted upload.

## Next precise target

The next support in exact support order is support `260`. Subsequent supports
in the current diagnostic ordering include `268` and `269`.

Full `O_13`, and hence full `SU(2)_13`, remains open.
