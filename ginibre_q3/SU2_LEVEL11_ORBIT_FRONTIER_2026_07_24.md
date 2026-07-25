# `SU(2)_11` and the resolved rank-six odd-orbit frontier

Date: 2026-07-24

## Status: resolved

Odd-level simple-current lifting gives

```text
full SU(2)_11 scalar GKS2*  <=>  scalar GKS2* in O_11.
```

The exhaustive exact residual census now proves the orbit-ring side. Hence
`SU(2)_11` satisfies scalar `GKS2*` for arbitrary word length, and every
partial-character coefficient is nonnegative.

The complete theorem and census are recorded in

```text
SU2_LEVEL11_FULL_GKS2_2026_07_24.md.
```

## Exact census

After support-overlap reduction and parity/residual decomposition, the complete
problem contains 27,962 residual regimes. Their exact partition is

```text
pointwise positivity                  2,882
direct Hall transport               22,443
braid transport                         371
Farey-fan transport                     248
translated tails and exact faces        241
AM-GM and cone certificates           1,700
additional AM-GM certificates            50
structural theorem sectors               27
                                      ------
total                                 27,962.
```

A literal key-set equality audit verified that these classes are disjoint and
exhaustive.

## Final obstruction

The last chamber was resolved by a uniform `0.14` log-base margin outside the
simplex of residual total below 75 and an exact five-thread fusion check of all
73,150 points inside that simplex. The minimum exact value was 2,738.

## Publication status

The exact full source-and-ledger archive was built and replayed but was lost in
a working-container reset before repository persistence. The theorem record
and PASS transcript are committed; restoring the one-command replay package is
a separate archival task. The previously published first-chamber verifier
remains reproducible directly from the repository.

## New odd-level frontier

The next odd-level target is

```text
O_13  <=>  full SU(2)_13.
```
