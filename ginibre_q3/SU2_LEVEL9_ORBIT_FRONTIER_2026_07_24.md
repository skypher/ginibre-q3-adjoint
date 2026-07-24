# `SU(2)_9` and the rank-five odd-orbit frontier

Date: 2026-07-24

## 1. Exact status

The odd-level simple-current lifting theorem proves

```text
full SU(2)_9 scalar GKS2*  <=>  scalar GKS2* in the rank-five orbit ring O_9.
```

Once the scalar theorem is known, adjoining one plus factor gives every
partial-character coefficient in the full fusion ring.  Thus there is no
separate full-ring chamber problem at level nine: the remaining finite
problem is exactly the orbit-ring problem.

Use the even-lift orbit basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6, B_4=V_8.
```

The earlier cyclic-gradient reduction identifies every orbit partial
coefficient with an adjacent Monge difference on `C_11 x C_11`.  Removing a
factor of largest radius reduces the scalar theorem to the structured cyclic
outer-Turan monotonicity target

```text
H_a(R) >= H_(a+1)(R),
```

with all radii in the remainder at most `a`.  This is the current conceptual
form of the level-nine problem.

## 2. Proved sectors

The separated-pole top ray is already proved for every odd rank by the
rational generating-function argument in Lemma 22H3S.  At rank five this
covers the chamber whose only active orbit labels are the two extreme radii,
both with minus sign and minimum odd exponent.  The bounded direct tadpole,
unfolded path, cyclic-moment, and rational-series calculations agree through
the recorded rectangle.

The exact rank-five bounded orbit search also checks all support-disjoint
signed words through twelve factors, including every partial-character
target.  Those calculations are evidence and implementation checks; the
separated-pole lemma is the unbounded theorem in its sector.

## 3. First obstruction after the separated-pole theorem

After declaring the separated-pole sector proved, the transport-fan analyzer
stops at

```text
support_code=59
parity_code=2
residual_code=5
variables=1,4.
```

Decoded in the orbit basis, the active signs are

```text
B_1 minus,  B_2 plus,  B_4 minus,
```

with minimum exponent parities even, odd, even, and with unbounded residual
powers in the two extreme active variables `B_1,B_4`.

The failure is only a failure of pairwise capacitated Hall transport after
all primitive rational two-dimensional fan rays of numerator and denominator
at most one hundred and every common shift through one hundred.  It is not a
negative `GKS2*` value.

This chamber is therefore the first exact all-exponent target beyond the
separated-pole theorem.  A proof may use a non-pairwise allocation, a Turan
identity coupling several spectral pairs, or a direct recurrence for its
integer coefficient sequence.

## 4. Negative searches already excluded

A deterministic tropical search tested one million signed exponent
directions at every orbit rank from five through fifteen, with direction
coordinates at most fifty.  It found no unique negative asymptotic maximizer.
This rules out a large bounded family of immediate high-power
counterexamples, but it is not a positivity proof because tied asymptotic
faces and untested directions remain.

## 5. Precise next theorem

The next finite theorem needed for the original program is:

> **Rank-five orbit theorem.** Every even-minus support-disjoint signed word
> in `O_9` has nonnegative scalar corner.

By odd simple-current lifting, this theorem is equivalent to full scalar
`GKS2*` in `SU(2)_9`; by the plus-factor identity it also yields the complete
partial-character column.  Combined with the sharp half-total stability
bound, it would prove every ordinary `SU(2)` partial coefficient with

```text
max(max_i a_i, c, ceil((T+c)/2)) <= 9.
```

The rank-five orbit theorem, rather than the former full level-five chamber
census, is now the first unresolved odd-level finite target.
