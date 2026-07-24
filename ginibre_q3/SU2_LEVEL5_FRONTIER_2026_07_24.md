# `SU(2)_5` finite-level frontier after the `k<=4` theorem

Date: 2026-07-24

## Status

The preceding note proves full all-length `GKS2*` and partial-character
positivity for `SU(2)_k`, `1<=k<=4`.  This note records the first complete
spectral chamber census for level five.  It is **not** a proof of the level-five
theorem.

Let `V_1,...,V_5` be the nontrivial simples.  After eliminating labels which
occur with both signs, each label is absent, plus, or minus.  Fixing the
parity of every positive exponent and requiring even total minus parity gives
exactly

```text
1441 support/parity chambers.                         (1)
```

Write each positive exponent as

```text
p_a=p_a^0+2r_a,
```

with `p_a^0` equal to one or two according to its parity.  The Verlinde
spectral formula then writes the corner as a finite signed sum of exponential
monomials in the residual variables `r_a`.

A high-precision chamber analyzer applies the same coordinatewise
capacitated-Hall transport used exactly at level four.  Its census is

```text
all chambers                 1441
proved by direct Hall        1312
not covered by direct Hall    129.                    (2)
```

Thus the elementary transport proves more than ninety-one percent of all
level-five parity chambers at arbitrary exponent depth.  The 129 failures
are failures of this sufficient transport criterion, not negative fusion-ring
coefficients.

Grouping the failures by their positive and negative spectral coefficient/base
arrays leaves

```text
67 distinct spectral signatures.                     (3)
```

Fifty signatures occur twice, thirteen occur once, and four occur four times.
Simple reversal of the five label coordinates alone is not the full symmetry:
it gives 115 orbits, so further compression must use the simple-current and
node symmetries at the spectral level rather than bare label reversal.

Seven non-Hall regimes already have only two active labels.  These provide a
small first target for exact closed-form analysis before attacking the
higher-dimensional signatures.

## Exact bounded fusion regression

The independent program

```text
character_ring_iter/verify_su2_level5_regression.py
```

evolves the integral `SU(2)_5` fusion ring directly.  With every active
exponent in `1,...,6`, it checks all support-disjoint sign chambers of even
minus parity and every target `V_c`, `0<=c<=5`.  The recorded result is

```text
SU2_LEVEL_5_REGRESSION PASS
cases=177243
partial_coefficients=1063458
minimum=0
minimum_positive=1
max_power=6.                                          (4)
```

This is exact bounded evidence for the full partial-character conjecture.  It
is not an all-exponent proof.

## Next exact targets

The efficient order of attack is now:

1. derive exact formulas for the seven two-label non-Hall regimes;
2. identify the full simple-current action on spectral signatures;
3. replace the numerical chamber census by exact algebraic-number comparisons;
4. search for a small collection of Turan or crossed-pair inequalities which
   covers the remaining signatures, as happened at level four.

No counterexample to `GKS2*` or partial-character positivity was found.
