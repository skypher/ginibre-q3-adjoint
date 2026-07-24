# `SU(2)_5` chamber census and its resolution

Date: 2026-07-24

## Current status: resolved

The original version of this note recorded a direct spectral chamber census
for the full level-five fusion ring and correctly marked 129 chambers as not
covered by coordinatewise Hall transport.  Those chambers are now resolved
by the odd-level simple-current lifting theorem in

```text
SU2_ODD_SIMPLE_CURRENT_LIFT_2026_07_24.md.
```

The rank-three simple-current orbit ring `O_5` already had an all-length proof
(Proposition 22H3N).  The exact lifting identity proves that scalar `GKS2*` in
`O_5` is equivalent to scalar `GKS2*` in the complete ring `SU(2)_5`.
Adjoining one plus factor then proves the full partial-character column.

Therefore:

```text
SU(2)_5 satisfies full all-length GKS2*,
and every partial-character coefficient is nonnegative.
```

## Historical direct-transport census

Let `V_1,...,V_5` be the nontrivial simples.  After eliminating labels which
occur with both signs, each label is absent, plus, or minus.  Fixing the
parity of every positive exponent and requiring even total minus parity gives

```text
1441 support/parity chambers.
```

Writing each positive exponent as `p_a=p_a^0+2r_a`, the Verlinde formula
turns each chamber into a finite signed sum of exponential monomials.  The
direct coordinatewise capacitated-Hall test gave

```text
all chambers                 1441
proved by direct Hall        1312
not covered by direct Hall    129.
```

The 129 failures collapse to 67 spectral signatures.  They are now known to
be failures only of that sufficient transport ansatz: the simple-current
lifting theorem proves their nonnegativity simultaneously.  In particular,
the earlier proposed program of deriving 67 separate spectral inequalities
is unnecessary for level five.

The direct census remains useful diagnostically because the same transport
ansatz reappears in higher odd orbit ranks.

## Independent bounded regression

The program

```text
character_ring_iter/verify_su2_level5_regression.py
```

checks the complete integral `SU(2)_5` fusion ring with every active exponent
in `1,...,6`, every support-disjoint even-minus chamber, and every target.
Its recorded output is

```text
SU2_LEVEL_5_REGRESSION PASS
cases=177243
partial_coefficients=1063458
minimum=0
minimum_positive=1
max_power=6.
```

This remains an implementation cross-check.  The all-length theorem is the
orbit proof plus exact simple-current lifting.

## New frontier

For odd levels, the full-ring problem is now exactly the orbit-ring problem.
The next unproved odd full fusion ring is therefore `SU(2)_9`, equivalently
the rank-five orbit ring `O_9`.  Existing files already prove its separated-
pole top ray and record the first transport obstructions beyond that sector.
