# A rank-parameterised AM-GM certificate stack for the odd orbit rings

Date: 2026-07-25

## Why

The bulk AM-GM stage was the largest single component of the lost `O_11`
proof: 1,558 certificates covering most of the hard residual regimes. It has
to be rebuilt, and the coverage scoring in
`SU2_O11_CENSUS_RECONSTRUCTION_2026_07_25.md` showed why sector theorems
cannot substitute for it. Of the 2,586 residual keys, only 14 are supported on
two labels; 2,433, or 94.1%, live on four or five labels.

Rather than rebuild it for rank six alone, this is built once and parameterised
by rank, so the same code serves `O_7`, `O_9`, `O_11`, `O_13` and beyond.

## 1. Certificate search is a linear program

For a regime, write the spectral sum as

```text
sum_t s_t c_t prod_l lambda_(t,l)^(r_l),        r_l >= 0,
```

with `s_t` a sign, `c_t>0`, `lambda_(t,l)>0`. Allocate a share `alpha[p][x]`
of each positive term `p` to each negative term `x`, and apply the weighted
AM-GM inequality using the allocation itself as the weight vector. The
requirement then splits into four families:

```text
(N)  sum_p alpha[p][x] = 1                              each negative x
(C)  sum_x alpha[p][x] <= 1                             each positive p
(D)  sum_p alpha[p][x] log lambda[p][l]
        >= log lambda[x][l]                             each x, free l
(K)  sum_p alpha[p][x] log c[p] >= log c[x]             each negative x
```

Every one is **linear in alpha**. So certificate search is a capacitated
transportation problem, not a general nonlinear search, and the "denominator
100 capacitated weighted AM-GM allocation" language in the earlier notes is
literally a rational point of that polytope.

The program maximises the worst-case slack, so a feasible answer arrives with
margin to spare for rational rounding.

## 2. Trust boundary

Floating point proposes; it never decides. A proposal is rounded to a fixed
rational denominator and then re-verified in 512-bit interval arithmetic, with
the spectral data rebuilt in interval mode from `2cos(pi t/n)` so that
verification cannot inherit the proposal's arithmetic. Structural conditions
`(N)` and `(C)` are checked exactly over the rationals. A certificate is
accepted only when every slack has a strictly positive interval lower bound.

Denominators escalate `100, 200, 500, 1000, 5000` and stop at the first that
verifies, so each certificate records its own cost.

## 3. Results

```text
rank  level  direct Hall  residual  closed  coverage
  4      7          191        13      10     76.9%
  5      9        2,207       201     169     84.1%
  6     11       22,491     2,589   1,998     77.2%
  7     13       26,377       400*    216     54.0%
```

`*` rank seven is a sample of the first 400 residual regimes in support
order, not a random draw, so its rate is indicative rather than an estimate
for the whole ring.  `O_13` has 34,924 failure keys in total.

Rank six completes in 258 seconds. Its denominator profile is

```text
100:1183  200:446  500:158  1000:122  5000:89.
```

The recorded `O_11` run used 1,558 bulk AM-GM certificates. This stack closes
**1,998** regimes, so the automated search is not merely reproducing the lost
stage but strictly exceeding it, cutting the set needing other machinery from
roughly 887 to 591.

Of the 591 remaining at rank six, 408 are genuinely LP-infeasible, meaning no
single allocation exists and the regime needs a decomposition, and 183 have a
feasible allocation whose margin did not survive rounding at any denominator on
the ladder.

If the sampled rate were to hold, the stack would close on the order of
nineteen thousand `O_13` regimes automatically.  The exact program currently
closes 211 by hand.  Even allowing for the sampling bias and for the rate
falling on harder supports, that is a change of scale rather than of degree,
and it is the argument for parameterising by rank rather than rebuilding per
level.

The rank seven sample took 636 seconds for 400 regimes, so the full ring is
roughly fifteen hours in Python.  That is tolerable once but is the point at
which a C++ port of the inner loop becomes worthwhile.

## 4. Validation

The Python enumeration reproduces the C++ `O_11` census exactly on all six
counts: 1,230 chambers, 25,076 regimes, 22,490 direct, 4 separated, 2,882
pointwise, 2,586 failures. At rank five it yields 2,408 regimes, matching the
`certified_regimes=2408` of the recorded `O_9` proof.

Soundness was checked independently of the certificate logic by evaluating the
actual signed sum at 7,200 exponent points across 120 certified regimes, with
no violation. As a negative control the verifier was fed allocations with
broken normalisation, broken capacity, and permuted columns, and rejected all
three.

Two real solver bugs were caught by making the LP validate its own output:
negating a row of `A x <= b` reverses that inequality rather than preserving
it, and an artificial variable left basic after phase one lets phase two pivot
around a row it cannot represent. Both produced proposals that violated their
own constraints, and both would otherwise have surfaced as mysterious
certificate failures.

## 5. Caveats

The residual count here is 2,589 against the census's 2,586. The direct Hall
test is a tolerance-based diagnostic, and three regimes fall on the other side
of it under mpmath rather than long double. This is the same class of
disagreement as the 47-regime discrepancy already recorded, and it affects
which regimes are *attempted*, never whether one is *accepted*.

The certificate found for a regime depends on the simplex pivot rule, since
degenerate vertices are common. Different rules yield different allocations of
equal margin, and occasionally one rounds successfully where another does not.
Acceptance is unaffected: only interval-verified allocations are counted.

This stack is a certificate engine, not a theorem. It closes 77.2% of the
rank-six residual regimes; the remainder still needs the decomposition
techniques, and `O_11` is not reproved until those are supplied.

## 6. Replay

```text
python3 character_ring_iter/verify_su2_orbit_amgm_stack.py --rank 6
```

Transcripts: `certificates/su2_orbit_amgm_rank4.log`, `...rank5.log`,
`...rank6.log`.
