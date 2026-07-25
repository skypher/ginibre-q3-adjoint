# `SU(2)_11` rank-six orbit frontier — resolved

Date: 2026-07-24

The former frontier is now closed. The rank-six simple-current orbit ring
`O_11` satisfies scalar `GKS2*` for arbitrary word length. By odd-level
simple-current lifting, full `SU(2)_11` satisfies scalar `GKS2*`, and every
partial-character coefficient is nonnegative.

The authoritative proof is

```text
SU2_LEVEL11_FULL_GKS2_2026_07_24.md
```

with deterministic replay through

```text
character_ring_iter/replay_su2_o11_full_gks2.sh.
```

The exhaustive exact census contains 27,962 residual regimes:

```text
pointwise       2,882
direct Hall    22,443
braid             371
post-braid      2,266
```

The post-braid archive is partitioned into 1,623 global AM-GM keys, 248
additional AM-GM keys, 387 translated-orthant trees, two uniform-margin finite
triangles, and six previously proved structural keys.

The ordinary stable consequence is now

```text
max(max_i a_i,c,ceil((T+c)/2)) <= 11.
```

The next odd-level frontier is `O_13`, equivalently `SU(2)_13`.
