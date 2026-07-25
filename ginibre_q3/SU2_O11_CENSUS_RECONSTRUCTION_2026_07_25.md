# Reconstructed `O_11` residual census

Date: 2026-07-25

## Why

`SU2_LEVEL11_FULL_GKS2_2026_07_24.md` records a complete exhaustive proof for
the rank-six orbit ring, partitioning 27,962 residual regimes across eight
certificate classes. The source-and-ledger package for that run was lost
before persistence, and both branch copies are corrupt: the archive on
`agent/o11-replay-restore-20260725-v2` is truncated at the writer's byte cap,
and `verify_su2_o11_amgm_frontier_exact.cpp` on `agent/o11-amgm-exact-20260724`
ends mid-literal in binary garbage.

So the level-eleven theorem currently has a transcript and no reproduction
path. Before any of it can be regenerated, the *census* has to exist again:
without it there is no way to say which regimes a new sector theorem covers,
or how much of the 27,962 remains.

This note reconstructs that census.

## Reconstruction

The rank-seven O13 analyzer is parameterised by rank, label count, and
Verlinde order. Rank six is the same program at

```text
kRank=6, kLabels=5, kOrder=13,
```

with the separated-pole signature on `B_1,B_5` and the support loop bounded by
`3^5=243`. `character_ring_iter/analyze_su2_o11_frontier.cpp` is that
specialisation, extended to count the pointwise chambers it would otherwise
skip.

## Agreement with the recorded run

```text
O11_DIRECT_CENSUS chambers=1230 regimes=25076 direct=22490
                  separated=4 pointwise=2882 failures=2586 threads=5
```

The analyzer reports only non-pointwise, non-separated regimes, so the
recorded total is recovered as

```text
25,076 + 2,882 + 4 = 27,962.                          (1)
```

Three independent quantities agree exactly with
`certificates/su2_o11_full_gks2_exact.log`:

```text
total      27,962   matches TOTAL
pointwise   2,882   matches POINTWISE
separated       4   matches the separated-pole sector.
```

## The 47-regime discrepancy

The fourth quantity does not agree exactly:

```text
direct Hall     recorded 22,443   reconstructed 22,490   (+47)
residual hard   recorded  2,633   reconstructed  2,586   (-47)
```

where the recorded hard count is `2,637` structural minus the `4` separated
regimes already counted. The two deviations are equal and opposite, so the
partition is internally consistent and the disagreement is entirely about
where 47 borderline regimes fall.

This is expected. The census is a long-double diagnostic with a relative
tolerance of `1e-12` in its Hall comparison, and both
`SU2_LEVEL13_ORBIT_FRONTIER_2026_07_25.md` and the level-eleven note state
that the direct census is a frontier diagnostic rather than proof evidence. A
regime that the diagnostic calls direct still has to be proved by an exact
certificate, so a tolerance that is slightly too generous costs coverage, not
soundness.

The 47 regimes are therefore a **worklist**, not a contradiction: they are
exactly the regimes where the diagnostic and the original exact classification
disagreed, and any regeneration must supply exact certificates for them
regardless of which side of the tolerance they fall on.

## Status

This is a diagnostic instrument, not a proof of anything. It does not restore
the level-eleven theorem. What it restores is the ability to measure: the
27,962-regime index against which the regeneration, and the newer sector
theorems on `B_1,B_5` and the top ray, can be scored.

Replay:

```text
g++ -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
    -Wsign-conversion -Wshadow -Werror -pthread \
    character_ring_iter/analyze_su2_o11_frontier.cpp -o analyze_su2_o11_frontier
./analyze_su2_o11_frontier
```

Transcript: `certificates/su2_o11_frontier_census.log`.
