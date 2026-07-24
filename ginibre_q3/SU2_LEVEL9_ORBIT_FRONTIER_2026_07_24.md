# `SU(2)_9` rank-five orbit frontier — resolved

Date: 2026-07-24

## Status

The rank-five simple-current orbit theorem is now proved. Consequently the complete fusion ring `SU(2)_9` satisfies scalar `GKS2*` for every signed word of arbitrary length, and every partial-character coefficient is nonnegative.

The proof and exact replay are in:

```text
SU2_LEVEL9_GKS2_2026_07_24.md
character_ring_iter/verify_su2_odd_orbit_rank5_exact.py
certificates/su2_odd_orbit_rank5_exact.log
```

The earlier transport obstruction

```text
support_code=59
parity_code=2
residual_code=5
variables=1,4
```

was not a negative coefficient. It was the first chamber not covered by pairwise Hall transport. The completed proof introduces endpoint-lattice certificates, which reserve positive capacity for the easy negative terms and cover the remaining hard term on two complementary integer half-planes.

## Exact census

The complete rank-five orbit calculation has:

```text
272 support/sign/parity chambers
 66 pointwise-nonnegative chambers
560 pointwise residual regimes
```

The other 2,408 residual regimes are exhausted by:

```text
direct coordinatewise Hall transport     1828
endpoint-lattice transport                 332
exact integral fusion leaves               205
shifted braid transport                      39
separated-pole theorem                        4
```

Every algebraic comparison is checked using rational intervals around Sturm-isolated roots of

```text
x^5+x^4-4x^3-3x^2+3x+1.
```

The certificate selector is not trusted for coverage or signs: the verifier reconstructs the complete chamber/key set and rechecks every selected inequality.

## Consequence

By the odd-level simple-current lifting theorem,

```text
O_9 scalar GKS2*  <=>  full SU(2)_9 scalar GKS2*.
```

Adjoining one plus factor gives the full partial-character column. Combining the finite theorem with sharp half-total stability proves every ordinary `SU(2)` partial coefficient satisfying

```text
max(max_i a_i,c,ceil((T+c)/2)) <= 9.
```

The next unresolved odd-level finite problem is the rank-six orbit ring `O_11`, equivalently the complete fusion ring `SU(2)_11`.
