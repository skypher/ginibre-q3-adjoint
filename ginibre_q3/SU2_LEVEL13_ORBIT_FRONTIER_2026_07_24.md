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

This census uses long-double spectral arithmetic and is a frontier diagnostic.
Every theorem below is instead replayed with rational Sturm isolation, exact
rational capacities, directed-rounding MPFR logarithms where needed, and
integral fusion arithmetic for finite leaves.

## Exact initial-frontier advance

The exact program now closes 172 direct-Hall failure keys.

### First two rays

The previously published first obstruction is

```text
support=56, parity=3, residual=1,
(B_1^-)^(1+2p) B_4^-,                 p>=0.
```

Exact leaves are `0,0,8`, and a denominator-100 weighted AM-GM certificate
proves the tail from residual floor two.

The next ray is

```text
support=83, parity=2, residual=1,
(B_1^-)^(2+2p) B_5^+,                 p>=0.
```

Two spectral pairs vanish identically. Exact values are `0,0,10`, and a
second denominator-100 AM-GM certificate proves the entire remaining tail.

### Complete `B_1^- B_5^-` chamber

All residual faces of

```text
(B_1^-)^(1+2p)(B_5^-)^(1+2q),         p,q>=0,
```

are proved. The two axes have exact initial values `0,0,2` and `0,0,20`,
respectively. The interior has a two-coordinate denominator-100 allocation
and first exact value eight. Thus the whole two-parameter chamber is
nonnegative.

### Decorated rays and support 169

The ray

```text
(B_1^-)^(1+2p) B_2^+ B_5^-
```

has exact leaf zero and a denominator-100 tail beginning with value two.
At support `169`, the ray `B_1^+B_2^-B_5^-` and both three-variable parity
interiors are also proved. Their first interior values are `7914` and `2348`.

### Generic multivariable replay

A reusable strict C++ verifier proves 163 additional three-, four-, and
five-variable residual regimes:

```text
support 178      8
support 187      8
support 196      8
support 205     32
support 214     32
support 223      8
support 232     32
support 241     32
support 250      3
                 ---
total           163.
```

The support-205 tables include six denominator-1000 certificates and two
near-boundary denominator-10000 tables. Directed rounding confirms every
geometric and capacity inequality. The complete 163-case replay runs in
about 4.84 seconds and uses under 6 MiB resident memory on the recorded host.

The exact key count is therefore

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

## Reproducibility package

The exact sources, allocation tables, and transcripts are stored in a
SHA-256-pinned compressed package split into GitHub-safe base64 chunks:

```text
certificates/su2_o13_initial_frontier_package.b64.part00
certificates/su2_o13_initial_frontier_package.b64.part01
certificates/su2_o13_initial_frontier_package.b64.part02
certificates/su2_o13_initial_frontier_package.b64.part03
```

The reconstructed archive must have hash

```text
a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2.
```

Replay with

```text
character_ring_iter/replay_su2_o13_initial_frontier_package.sh
```

The run transcript is

```text
certificates/su2_o13_initial_frontier_exact.log.
```

## Next precise target

The next unresolved key in support order is

```text
support=245, parity=2, residual=1.
```

It decodes to the one-variable tail

```text
(B_1^-)^(4+2r) B_6^+,                  r>=0.
```

The other two non-Hall residual faces at the same support have residual codes
`2` and `3`. This support is the next exact rank-seven orbit target.

Full `O_13`, and hence full `SU(2)_13`, remains open.
