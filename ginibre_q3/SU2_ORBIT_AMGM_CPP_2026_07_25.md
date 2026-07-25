# The AM-GM stack in C++, and the first bulk `O_13` result

Date: 2026-07-25

## Why port

The Python stack in `SU2_ORBIT_AMGM_STACK_2026_07_25.md` is correct but slow.
A full rank-six pass took 258 seconds and a full rank-seven pass was projected
at roughly fifteen hours, so every measurement cost more than the increment it
was measuring. Measurement, not mathematics, had become the bottleneck.

## Result

Full `O_13`, the rank-seven orbit ring, in **41 seconds**:

```text
SU2_ORBIT_AMGM_CPP rank=7 level=13 threads=5
  direct_hall      208085
  residual          34927
  pointwise         14896
  lp_feasible       29238
  certified         22596
  denominators      100:10728 200:5688 500:2837 1000:2015 5000:1328
  coverage          64.7% of residual regimes
```

The exact hand-built program currently closes **211** `O_13` failure keys. This
closes **22,596** of 34,927, in well under a minute.

All ranks:

```text
rank  level  residual  certified  coverage   wall
  4      7        13         11     84.6%   <1s
  5      9       201        185     92.0%   <1s
  6     11     2,589      2,052     79.3%    2s
  7     13    34,927     22,596     64.7%   41s
```

Rank six is about 250 times faster than the Python at the same task, and the
C++ certifies slightly more at every rank, 2,052 against 1,998 at rank six,
because long double logarithms give the linear program better conditioned data
than Python's float conversion did.

This is the flat allocation stage only. The recursive orthant decomposition
that lifted Python's rank six to 87.4% is not yet ported, so these numbers are
a floor rather than the technique's ceiling.

## Agreement with the recorded census

The enumeration reproduces the published `O_13` diagnostic census:

```text
                 recorded    C++
direct Hall       208,084   208,085
non-Hall           34,924    34,927
```

Off by one and by three respectively, the same tolerance-class disagreement
already recorded for `O_11`, where the long double Hall comparison places a few
borderline regimes on the other side of its threshold. It changes which regimes
are attempted, never which are accepted. At rank five the pointwise count is
560, matching the recorded `O_9` `pointwise_regimes=560` exactly.

## Trust boundary

Unchanged from the Python. Long double proposes an allocation by linear
programming and nothing it produces is accepted. Every certificate is rounded
to a rational denominator and re-checked in 512-bit MPFR interval arithmetic
with directed rounding, on spectral data rebuilt independently in interval
mode; the normalisation and capacity conditions are checked exactly over the
integers.

Two spectral models are therefore built from scratch by different arithmetic,
and the program cross-checks them at startup, requiring agreement to long
double precision at every node and label. That check must compare at long
double precision rather than demand containment: the interval enclosures are
tight to about `1e-154`, while the long double values carry their own error of
about `1e-19`, so a float value legitimately lies outside its own enclosure.

## Validation

Every recorded run carries its own soundness sample and negative control:

```text
rank  soundness points  violations   corrupted  refused
  4              440             0          22       22
  5            7,400             0         370      370
  6           71,800             0       1,000    1,000
  7           80,000             0       1,000    1,000
```

Soundness evaluates the actual signed sum at pseudo-random exponent points in
each sampled certified regime, independently of the certificate logic. The
negative control corrupts a verified allocation, breaking normalisation and
then capacity, and requires the verifier to reject it.

The soundness check earned its place immediately: its first run reported 690
violations in 12,000 points. The cause was the sampler, not the certificates --
it advanced its random state inside the term loop, so every term was evaluated
at a different exponent vector rather than all terms at one shared point. With
one point per trial the violation count is zero at every rank. A checker that
had merely been assumed correct would have reported nothing either way.

## Caveats

Coverage falls as rank rises, 92.0% at rank five to 64.7% at rank seven. Higher
rank means more spectral terms and more free coordinates, so a single
allocation covers proportionally less. That is the expected shape and it is the
argument for porting the decomposition next, not evidence against the method.

The 12,331 rank-seven regimes still open are the real remaining work. They are
now enumerated and addressable in seconds rather than hours, which is the point
of this increment.

## Replay

```text
make -f character_ring_iter/Makefile.research verify_su2_orbit_amgm_stack
./character_ring_iter/verify_su2_orbit_amgm_stack --rank 7 --threads 5 \
    --soundness 400 --control
```

Transcripts: `certificates/su2_orbit_amgm_cpp_rank4.log` through `rank7.log`.
