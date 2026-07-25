# What the surviving orbit regimes actually need

Date: 2026-07-25

## Why measure before building

The order-cone fan was built on the reasonable guess that the surviving
regimes were diagonally obstructed. It closed exactly one regime at rank six.
That is a cheap lesson but a real one, so this note characterises the open
regimes before any further technique is written.

## Method

`verify_su2_orbit_amgm_stack --dump-open` emits one line per regime that
survives both the flat allocation stage and the recursive orthant
decomposition, recording its free-coordinate count, its negative and positive
term counts, the split nodes spent on it, and, crucially, whether a root
allocation existed at all.

That last field separates two failures that look identical in a coverage
number and need opposite remedies:

* **root LP infeasible** -- no capacitated allocation exists. The regime needs
  a genuinely different technique.
* **root LP feasible** -- an allocation exists and the linear program found it,
  but the certificate did not verify after rational rounding. The mathematics
  is fine; the rounding lost the margin.

## Result at rank six, 326 open regimes

```text
root LP status
  infeasible   168   51.5%
  feasible     158   48.5%

free coordinates
  k=1     25    7.7%
  k=2     56   17.2%
  k=3    115   35.3%
  k=4    100   30.7%
  k=5     30    9.2%

negative terms
  6:  37     7: 161     8:  88     9:  40

k against root LP
  k=1  infeasible 25   feasible  0
  k=2  infeasible 43   feasible 13
  k=3  infeasible 45   feasible 70
  k=4  infeasible 41   feasible 59
  k=5  infeasible 14   feasible 16
```

50 distinct supports, the largest contributing 24 regimes, so the failures are
spread across the chamber structure rather than concentrated in a corner of it.

## What this says

**Nearly half the open regimes do not need a new technique.** For 158 of 326 an
allocation exists and was found; the certificate was lost in rounding. That is
a defect in how the proposal is converted to a rational, not in the method, and
it should be recoverable by making the linear program aware of the rounding it
will subsequently undergo -- requiring a margin large enough to survive
perturbation of each weight by up to `1/2d`, rather than maximising the slack
and hoping. The denominator ladder already climbs to 5,000 and does not rescue
these, which is consistent with a vertex that is thin in a direction the ladder
cannot fix.

**Every one-coordinate survivor is genuinely infeasible**, 25 of 25. A ray with
no allocation cannot be helped by finer rounding, and the orthant recursion
already walks it to depth 240. These need a different argument.

**The mass sits at three and four free coordinates**, 66% together, not at the
extremes. The fan was tested against this population and closed one regime,
which is now explained: a fan only helps when the obstruction is diagonal, and
these are mostly not.

## Consequence for the plan

The next increment is rounding-aware margin selection, targeting the 158
feasible-but-unverified regimes at rank six and the corresponding population at
rank seven. It is a small change to the linear program rather than a new
certificate class, and the measurement says it addresses about half the
remaining work.

The 168 infeasible regimes at rank six remain genuinely open. The node budget
is never the binding constraint at any rank measured, so they will not yield to
more compute.

## Replay

```text
./character_ring_iter/verify_su2_orbit_amgm_stack --rank 6 --decompose --dump-open
```

Transcript: `certificates/su2_orbit_open_rank6.log`.
