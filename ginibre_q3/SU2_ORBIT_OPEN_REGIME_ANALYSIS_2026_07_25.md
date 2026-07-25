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

## Follow-up: the diagnosis was wrong, the fix worked anyway

The analysis above predicted that the root-feasible regimes were losing their
certificates to rational rounding, and that a rounding-aware margin would
recover them.  Measuring the linear program's worst-case slack on exactly those
regimes refutes it:

```text
root margin, root-feasible open regimes at rank six
  1e-4 .. 2e-3      2    1.3%
  2e-3 .. 1e-2      4    2.6%
  >= 1e-2         145   96.0%
```

The perturbation scale at the finest denominator, 5,000, is `2e-4`.  Not one
regime has a margin below it, and 96% sit two orders of magnitude above.  The
margins were never thin, so "lost in rounding for want of resolution" was
simply the wrong explanation.

Replacing round-to-nearest-plus-fixup with a capacity-aware greedy -- start
from the floor and hand out each remaining unit to whichever positive term
maximises the resulting minimum slack -- nonetheless produced a large gain:

```text
rank 7, flat allocation stage only
  round-to-nearest   22,596   64.7%   41s
  greedy             29,175   83.5%   39s
```

The reason is not resolution but structure.  Round-to-nearest is blind to the
constraints: a comfortable margin at the linear program's vertex says nothing
about whether the nearest lattice point of the allocation polytope is feasible,
because the slack moves by `sum_p |delta_p| |log lambda_p|` and the logarithms
are large.  The greedy optimises the quantity that actually decides the
certificate.

At rank six the same change moves 180 regimes from needing a split to closing
flat, 2,052 to 2,232, while the total rises only from 2,263 to 2,270 because
decomposition was already catching most of them.  The value is therefore mostly
in cost rather than coverage at rank six, and in both at rank seven, where the
greedy flat pass alone beats the previous flat-plus-decomposition result by 592
regimes while running twenty-four times faster.

## What blocks the genuinely infeasible regimes

For a regime whose root allocation is infeasible, re-solving with each
constraint family dropped in turn identifies the obstruction:

```text
rank six, 168 root-infeasible regimes
  domination (D) blocks                 89   53.0%
  both block individually               65   38.7%
  constant condition (K) blocks         14    8.3%
```

The one-coordinate survivors invert this. Of the 25 at rank six, thirteen are
blocked by the constant condition and only two by domination. On a single ray
there is one growth direction and it is usually dominable; what fails is paying
the negative term's magnitude.

## The indicated generalisation

The certificate as formulated requires each negative to be covered by a
*convex* combination of positives, `sum_p alpha[p][x] = 1`. That is stronger
than necessary. Weighted AM-GM needs only that the weights used in the
geometric mean sum to one, so allocate total mass `S >= 1` and set
`beta = alpha / S`:

```text
sum_p alpha_p v_p = S sum_p beta_p v_p >= S prod_p v_p^(beta_p).
```

The conditions become

```text
(D)  sum_p beta[p][x] log lambda[p][l] >= log lambda[x][l]
(K)  log S + sum_p beta[p][x] log c[p]  >= log c[x]
```

with capacity `sum_x S_x beta[p][x] <= 1`. The domination condition is
unchanged, and the constant condition gains a free `log S`. That is exactly the
slack the K-blocked regimes are missing, and it strictly subsumes the present
certificate, which is the case `S = 1`.

For fixed `S` the system is still linear, so the existing solver handles it by
sweeping a few values of `S` rather than by any new machinery. This targets the
14 K-blocked regimes and the 13 K-blocked one-coordinate rays directly, and may
reach some of the 65 that block on both.

It will not help the 89 blocked purely on domination. Those need the positives
to grow at least as fast as the negative in every free coordinate, which no
reweighting can manufacture.

## Replay

```text
./character_ring_iter/verify_su2_orbit_amgm_stack --rank 6 --decompose --dump-open
```

Transcript: `certificates/su2_orbit_open_rank6.log`.
