# `SU(2)_11` and the rank-six odd-orbit frontier

Date: 2026-07-24

## Exact reduction

By odd-level simple-current lifting,

```text
full SU(2)_11 scalar GKS2*  <=>  scalar GKS2* in the rank-six orbit ring O_11.
```

Once the scalar theorem is proved, adjoining one plus factor gives the complete partial-character column. Thus the next unresolved odd-level problem is entirely the rank-six orbit ring.

Use the even-lift basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6, B_4=V_8, B_5=V_10.
```

The cyclic-gradient and outer-Turan reductions from the odd-orbit program apply unchanged, now on `C_13 x C_13`.

## First all-exponent transport obstruction

After support-overlap reduction and removal of the already proved separated-pole sector, the generalized transport-fan analyzer stops first at

```text
support_code=83
parity_code=2
residual_code=3
variables=1,5.
```

Decoded in the orbit basis, this is the two-label chamber

```text
(B_1^- )^(2+2p) (B_5^+ )^(1+2q),      p,q>=0.        (1)
```

Thus the first rank-six obstruction is simpler in support than the first rank-five obstruction: only the two extreme active orbit labels occur. It is not covered by direct coordinatewise Hall transport or by the endpoint-reserve template that completed the rank-five theorem.

This is a failure of those sufficient transport ansatzes, not a negative coefficient.

## Exact bounded scan

The program

```text
character_ring_iter/verify_su2_o11_first_chamber.py
```

uses exact integer orbit-fusion dynamic programming. It evolves the two commuting paired transfer operators

```text
(B_1 tensor 1-1 tensor B_1)^2,
(B_5 tensor 1+1 tensor B_5)^2
```

from the base word `(B_1^-)^2 B_5^+`. The recorded rectangle is

```text
0<=p<=100,       0<=q<=100.
```

Its result is

```text
SU2_O11_FIRST_CHAMBER PASS
checks=10201
minimum=0
zeros=[(0,0),(0,1),(0,2),(1,0)].
```

Every other coefficient in the rectangle is strictly positive. This is bounded evidence only; the all-exponent two-variable inequality remains open.

The transcript is `certificates/su2_o11_first_chamber.log`.

## Precise next target

The next useful lemma is:

> For all nonnegative integers `p,q`, the scalar corner of (1) is nonnegative.

Because the two paired transfer matrices commute on the 36-dimensional tensor-square orbit space, the bivariate generating function is rational. A promising route is therefore to compute its exact denominator and prove coefficient positivity by a finite recurrence or a positive partial-fraction grouping. The four boundary zeros above must be built into any such formula.

A proof of this first chamber will allow the transport-fan census to advance to the next genuinely new rank-six obstruction. A proof of all rank-six chambers would, by simple-current lifting, establish full `GKS2*` and every partial-character coefficient for `SU(2)_11`, extending the ordinary stable range to `K(c)<=11`.
