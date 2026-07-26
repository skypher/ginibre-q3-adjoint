# Engine-hard families: stable zeros and candidate wall families

Date: 2026-07-26

## Measurement

The engine records which certificate class closes each regime
(`--dump-certified`).  The heavy classes — Farey merges, Abel chains, exact
integer leaves — mark regimes that are hard for this certificate engine.
They are candidates for the actual equality and wall structure, but the
certificate class alone does not prove mathematical tightness.  Aggregating
them by signed support family across the measured levels gives:

```text
level   hard regimes   families   containing V_k
  3            4            2            2
  4           16            8            7
  5           61           25           22
```

The natural guess — every tight family touches the affine wall — is
**false**, and its counterexamples are the finding:

```text
level 4:  V1- V3-                       (max label k-1... in fact 3 = k-1)
level 5:  V1- V3-   the same family     (max label 3 = k-2)
level 5:  V2- V4-   its double
level 5:  V1- V4+
```

`V1- V3-` is tight at level 4 and again, unchanged, at level 5.

## Why: these are the stable ordinary zeros

For a support-disjoint two-minus pair `a != b`, the subset formula gives the
minimal corner

```text
J = 2 N(a,b) = 0,
```

an exact zero **in the ordinary ring**. By the finite-to-stable theorem the
finite coefficient equals the ordinary one once `k >= a+b`, so this zero is
reproduced identically at every level from `a+b` on. The recurring tight
families with small labels are precisely these stable zeros; they are
level-independent by the stability theorem, not by accident.

## The measured pattern

The data suggest the candidate decomposition

```text
tight(k) = { stable ordinary zeros, present for all k >= their threshold }
           union
           { wall families containing V_k, specific to the level }.
```

This is verified only for the displayed engine-hard families at levels
three through five: all small-label families there are stable zeros, and
every remaining family contains the top label (22 of 25 at level 5, 7 of 8
at level 4, 2 of 2 at level 3).  It is not yet a theorem for arbitrary
levels or arbitrary equality families.

## Consequence for the factor-axis induction

This splits the induction's burden in two:

1. **Ordinary equality cases**: the stable zeros are visible entirely in the
   ordinary ring, where the subset formula makes them enumerable per word
   shape. Any all-`N` theorem must hold with equality exactly there, which
   constrains the form of a valid inductive invariant — it must be sharp on
   `J = 2N(a,b) = 0` pairs and their decorations.
2. **Wall corrections**: everything else concentrates at `V_k`, which is the
   pair-reservoir domain the six-factor proof already handles at the wall.
   The induction's genuinely new content is only there.

The corrected level-six certificate is complete, but its full
`--dump-certified` aggregation and the level-eight measurement have not yet
been recorded.  Those measurements test whether the candidate decomposition
persists beyond the levels where it was first observed.

## Replay

```text
./verify_su2_orbit_amgm_stack --level 5 --decompose --dump-certified \
    | python3 character_ring_iter/mine_su2_tight_cases.py 5
```
