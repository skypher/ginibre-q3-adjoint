# Ginibre Q3 author self-audits

These files were created and are maintained inside the author's repository.
They are self-audits and algorithmic redundancy controls, not independent
peer review.  The historical directory name is retained only to preserve
paths recorded in manifests and transcripts.  These files do not alter any
replay-pinned source, replace the clean 71-task fleet collection gate, or
serve as premises of a theorem.

`GINIBRE_1970_SCOPE_SOURCE_AUDIT.md` records the direct primary-source basis
for the phrase “full Q3”: all finite-family/sign quantifiers and Proposition
1 cone promotion are matched, while Example 4's full positive-definite cone
is a distinct and false nonabelian claim.

`HKKO_BOUNDED_LITTLEWOOD_SOURCE_AUDIT.md` records a direct page-level check of
the primary bounded-Littlewood source.  It verifies the three specializations
used for types B, C, and D, including every determinant sign, index shift, and
factor of two.

`RAINS_CLASSICAL_SCHUR_SOURCE_AUDIT.md` directly checks the orthogonal and
symplectic Schur-integral criteria and derives the even-orthogonal
determinant-twisted correction used immediately upstream of those
determinants.

## Exact moment and EGF formula

`verify_full_q3_moment_formula_independent.py` formally multiplies
`(X-Y)^(2a)(X+Y)^n` and checks every coefficient against the manuscript's
double moment sum.  It then compares that sum and the derivative exponential
generating function against direct exact expectations for four rational
finite laws, checks three arbitrary formal moment ledgers, and verifies the
index/sign skeleton in all eight load-bearing C++ implementations.  The
manuscript's boundedness and absolute-convergence justification is required,
and five adversarial mutations must be rejected.  Every calculation uses
Python integers and `Fraction`; no floating-point arithmetic is used.

From the repository root, run:

```sh
python3 \
  ginibre_q3/referee_audits/verify_full_q3_moment_formula_independent.py
```

## Separately implemented classification coverage self-audit

`verify_full_q3_classification_coverage_independent.py` parses the Part III
theorem statements and the independent 70-row low-classical ledger.  It
reconstructs the complete type-A, classical, and exceptional rank/odd-degree
case split; recomputes every stated residual-box cardinality; distinguishes
the harmless middle/high overlap at `D296` from a gap; and checks the final
finite-tuple/sign quantifiers and the all-positive-definite scope boundary.
Five in-memory adversarial mutations must be rejected.  This audit is
structural: it complements, but does not replace, the exact arithmetic
replays.

From the repository root, run:

```sh
python3 \
  ginibre_q3/referee_audits/verify_full_q3_classification_coverage_independent.py
```

## Parts I--II two-minus interface

`verify_full_q3_two_minus_interface_independent.py` checks the only imported
theorem row in Part III.  It parses both displayed Haar integrals, verifies
that substituting `a=1` gives the Parts I--II functional exactly, and compares
the compact-connected/simple group scope, `n >= 0` range, Haar
normalizations, central-cover equalities, theorem number, bibliography entry,
and replay/Part III manifest bindings.  Five in-memory interface mutations
must be rejected.  This is an interface audit, not a substitute for the
200-stage arithmetic replay proving the imported theorem.  Its optional
`--root` argument permits the same audit of the executed baseline snapshot and
the exact post-rebind publication tree.  The live publication source is an
explicit authorized state; prose or metadata changes must update and repin
this checker rather than leaving its source-identity test stale.

From the repository root, run:

```sh
python3 \
  ginibre_q3/referee_audits/verify_full_q3_two_minus_interface_independent.py
```

## Full cone promotion algebra

`verify_full_q3_cone_promotion_independent.py` audits the step from the single
adjoint generator to arbitrary finite tuples from its closed multiplicative
convex cone.  In the formal variables `u=X+Y` and `v=X-Y`, it reconstructs
every signed monomial with nonnegative binomial coefficients and verifies that
the parity of the total `v` exponent is exactly the parity of the number of
minus signs.  It also checks constants, both signs of Ginibre's product
factorization, all final theorem quantifiers, and the uniform-closure passage.
Five adversarial mutations must be rejected.

From the repository root, run:

```sh
python3 \
  ginibre_q3/referee_audits/verify_full_q3_cone_promotion_independent.py
```

## Exact SO(3) scope obstruction

`verify_full_q3_so3_obstruction_independent.py` reconstructs the counterexample
that rules out replacing the adjoint-generated cone by all real continuous
positive-definite functions.  It derives the degree-two and degree-four
invariant dimensions from the spin-one Clebsch--Gordan recurrence, enumerates
the three pairing tensors, inverts their exact Gram matrix, and independently
expands all 16 terms of the `--++` product.  It also checks the manuscript's
matrix-coefficient positive-definiteness identity and rejects five mutations.

From the repository root, run:

```sh
python3 \
  ginibre_q3/referee_audits/verify_full_q3_so3_obstruction_independent.py
```

## Full-scale H25--H28 B-frontier determinant oracle

`verify_frontier_bounded_littlewood_oracle.cpp` checks the 120 accepted
H25--H28 B-frontier bad-count totals by a global algorithm different from the
fleet programs' reverse-Pieri state traversal.  It reconstructs the finite
rank moment by the bounded-Littlewood determinant formula and checks the exact
identity

```text
bad(q,j) = stable(j) - moment(B_(q-1),j).
```

The archived wrapper deliberately included the then-current, separately
audited Part III implementation
`verify_full_q3_bcd_bounded_littlewood_gmp.cpp`; consequently this is
algorithmically independent of the fleet traversal, but it is not a second
source implementation of the determinant formula.  Those exact included
bytes are preserved under the distributed execution-source snapshot and are
the bytes pinned by the oracle manifest.  The live implementation has since
expanded its row ledger, while retaining the same determinant engine.  The
wrapper supplies its own B14--B46 ledger, stable-moment reconstruction bound,
accepted-log parser, and B-only evaluation path.

The archived run used machine C (`nb1cb2f`), an AMD Ryzen 9 7950X with 32
logical CPUs, g++ 13.3.0, and 24 OpenMP threads.  It ran at niceness 10 from
approximately `2026-07-16T00:51:56Z` through `2026-07-16T00:56:20Z`.  The
strict compilation flags were

```text
-O3 -march=native -std=c++20 -fopenmp
-Wall -Wextra -Wpedantic -Wconversion -Werror
```

with GMP linked by `-lgmpxx -lgmp`.  The machine-C binary SHA-256 was
`3431f3dbdd335096e3f026990523e0aab094b590186b4fdb1a0904178aa482fe`.
The exact wrapper also completed a full degree-121
AddressSanitizer/UndefinedBehaviorSanitizer run on machine C.  That build used
O2, 16 OpenMP threads at niceness 15, and
`ASAN_OPTIONS=detect_leaks=1:halt_on_error=1`.  All 22 primes and all 120
comparisons passed in 1,931 seconds with no sanitizer diagnostic.  Its binary
SHA-256 is `e13f9ba8...05ca2`, and its transcript SHA-256 is
`9f4f101f...22fd8`.  GCC 13 emitted `-Wmaybe-uninitialized` only from the
system libstdc++ `<regex>` implementation under instrumentation; that single
system-header diagnostic was demoted from `-Werror`.  The project source
retained every strict warning flag, and its ordinary O2/O3 builds pass with
all warnings fatal.  The exact flags and run identities are recorded in
`frontier_bounded_littlewood_oracle_sanitizer_machine_c.txt`.

The run used 22 distinct prime fields.  Their 682-bit CRT modulus strictly
exceeds the 663-bit maximum stable-moment bound, so every nonnegative finite
moment has a unique exact reconstruction.  All 120 bad-count values matched;
the transcript ends with `failures=0`, `ALL PASS`, and `__EXIT_STATUS=0`.
`verify_frontier_bounded_littlewood_transcript.py` independently enforces the
120-row scope, all exact subtractions, prime-index coverage, the CRT-bound bit
lengths, and unique terminal records.  The ordinary and sanitizer transcripts
have identical ordered 120-row arithmetic projections.  Five adversarial
tests are rejected: a changed accepted value, an omitted row, a duplicate
prime, a changed CRT modulus, and an unknown record.

From the repository root, the portable reproduction commands are:

```sh
g++ -O3 -march=native -std=c++20 -fopenmp \
  -Wall -Wextra -Wpedantic -Wconversion -Werror \
  ginibre_q3/referee_audits/verify_frontier_bounded_littlewood_oracle.cpp \
  -lgmpxx -lgmp -o /tmp/verify_frontier_bounded_littlewood_oracle
OMP_NUM_THREADS=24 nice -n 10 \
  /tmp/verify_frontier_bounded_littlewood_oracle \
  ginibre_q3/certificates/post_m29 \
  > /tmp/frontier_bounded_littlewood_oracle.log
python3 \
  ginibre_q3/referee_audits/verify_frontier_bounded_littlewood_transcript.py \
  /tmp/frontier_bounded_littlewood_oracle.log
```

The accepted inputs, sources, archived output, and transcript checker are
bound by `frontier_bounded_littlewood_oracle_manifest.sha256`.
