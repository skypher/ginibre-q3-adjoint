# Full adjoint-generated Q3 certificates

This directory contains exact and directed-rounding certificates used only by
`full_q3_extension.tex`.  They are separate from the accepted exactly-two-minus
certificate boundary of `paper.tex`.

`fullq3complete0002_optimus.log` is the hardened aggregate current-source
replay of the complete Part III target `make full-q3-extension`.  It ran on
`optimus` from `2026-07-15T12:45:11Z` to `2026-07-15T13:30:06Z` with two
low-priority OpenMP workers.  The transcript authenticates the final acyclic
Part III document and the executed prepublication source snapshot

```text
fb54ba33c1711ba404f368ce3af8def226d90facf56953958e2791068178b210  Makefile
48d269b8d0ab67d15e49f99256fd0e36dd559691bf2039d2c373ca99fc483c02  full_q3_extension.tex
a92f75a0cfdd5ea7f7bdea8eafc11a3481ff4227d2f29b7de6f156bcc1cdc26e  verify_full_q3_document.py
e7a46ec33608a9c54aca3b6154b4cd6b864d2f9b2ffc00117384edd331e9ac61  run_fullq3complete0002_optimus.sh
c147434141f42484b6a46364803f784e1c15a8d48d06954ce903938360f52f87  full_q3_source_manifest.sha256
```

The initial publication-only rebind in `paper.tex` produced source-manifest
SHA-256
`1cbd33b571e990c8b8bce47a1bd549129dea6cc89d965ec5951ec6f8f265fb59`.
The live final-source manifest, including the compressed modular checker, has
SHA-256
`cb7f417576f6bd94686883abbff70f8f81e797ed501ae6058771d321c19e52a2`.
The historical aggregate transcript and the archived execution snapshot keep
the prepublication hash above; the live rebind verifier authenticates mutable
publication files through the current manifest while requiring every archived
arithmetic source to remain byte-identical to the execution snapshot.

The same transcript records successful acyclic document, artifact,
stable-hierarchy, type-A,
exceptional, classical analytic, exact residual, directed-tail, and
bounded-Littlewood stages.  In particular it reproduces the exact scopes
`total_pairs=4869`, `rows_checked=58`, and `pairs=12993`; rebuilds the
30-page PDF twice; and ends with `__EXIT_STATUS=0`.  The transcript SHA-256 is

```text
47321d1bcd291d42e1c05473fc042de0243994743888070336d3983d9060b430.
```

Its 4,377-line arithmetic projection agrees with
`fullq3complete0001_machine_c.log` after normalizing only the recorded OpenMP
thread count from four to two.  The older transcript remains the independent
machine-C audit trail.  The adjacent SHA-256 sidecar authenticates the new
imported bytes.  The individual source-matched certificates below remain the
load-bearing arithmetic records;
the aggregate transcript proves that the final Make target consumes them and
recomputes every active Part III verifier from the pinned source snapshot.

`fullq3bcd0002_machine_c.log` is the source-matched raw C++/GMP replay for the finite residual
boxes in `B_18..B_21` and `D_31..D_70`.  It was compiled and run on machine C
(`nb1cb2f`, 16 OpenMP workers) from the source hashes recorded in
`FULL_Q3_BCD_RESIDUAL_RUN.md`.  Acceptance requires all of the following:

- the shared stable-moment recurrence check succeeds;
- every required reverse-Pieri moment correction is regenerated;
- the exact scope is `B_pairs=137`, `D_pairs=4732`, `total_pairs=4869`;
- every hierarchy value is positive; and
- the terminal marker is
  `FULL_Q3_B18_D31_RESIDUAL VERIFICATION: ALL PASS`.

The adjacent `.sha256` file authenticates the raw transcript.  The replayable
source, rather than this stored output, remains the mathematical certificate.

`fullq3bcd0001_machine_c.log` is the first successful run.  It has the same
scope and result but predates the standalone verifier refactor and lacks an
explicit wrapper exit marker, so it is retained only as an audit trail and is
not the cited certificate.

`fullq3bcdlowtail0001_optimus.log` is the complete directed-rounding replay for
the 58 remaining exact-MGF tail rows.  It was compiled with strict warnings
and run on `optimus` from these sources:

```text
a489eaf1ae72bd1207b643dba2bab75121e6a46a16e6de01dd428f33681f1014  verify_full_q3_bcd_low_tail_mpfr.cpp
643f31317ba1be2014db2a9334654499c38981c5a753fd948db177cb4093cfb8  post_m29_bc_layered_mgf_mpfr.cpp
7e913141d86e3f52d003b8cf5a1f2a0582d4a93722ecc60c04201b062acdabb9  full_q3_bcd_remaining_data.hpp
```

The transcript records 384-bit outward MPFR rounding, eight recovery layers,
180 Bessel-series recurrence updates plus a geometric remainder, the full
exact rational parameter schedule for every row, `rows_checked=58`, and the
terminal marker `FULL_Q3_BCD_LOW_TAIL VERIFICATION: ALL PASS`.  Its adjacent
SHA-256 file authenticates the raw output.

`fullq3bcd0010_machine_c.log` is the accepted complete exact
bounded-Littlewood replay
for all residual boxes left in `B_2..B_17`, `C_2..C_28`, and `D_4..D_30`.
It was compiled with strict warnings and run on machine C (`nb1cb2f`) with 24
OpenMP workers from these sources:

```text
257bcd0deaffea437ad25d2c76422b465ad33d7c2669df3d0d3b72c84134223d  verify_full_q3_bcd_bounded_littlewood_gmp.cpp
7e913141d86e3f52d003b8cf5a1f2a0582d4a93722ecc60c04201b062acdabb9  full_q3_bcd_remaining_data.hpp
```

The determinant series are evaluated in exact prime fields and reconstructed
by GMP CRT only after the modulus product exceeds the elementary character
bound.  The transcript checks stable prefixes, the `B_2/C_2` isomorphism,
moment bounds, positivity of exactly 12,993 hierarchy values, and the terminal
marker `FULL_Q3_BCD_BOUNDED_LITTLEWOOD VERIFICATION: ALL PASS`; the wrapper
then records exit status zero.  The same transcript also reruns all six exact
polynomial tails and all six exact rational-cap tails, with positive exact
GMP margins and an asserted `6 + 6` scope.  Its adjacent SHA-256 file
authenticates the raw transcript.

The final source also validates the consecutive 70-row rank sequence, all
coverage endpoints, the method counts `6 + 6 + 58`, and each rational-cap
onset before arithmetic begins.  Its 70 printed row counts independently sum
to 12,993.  `fullq3bcd0009_machine_c.log` is arithmetic-identical in every
moment and exact-tail margin but predates these fail-closed row-ledger checks;
it is retained as the immediately preceding audit trail.

`fullq3bcd0008_machine_c.log` is the immediately preceding residual-only
replay.  It has the same 12,993 residual result but predates integration of
the twelve exact-tail ledgers, so it is retained as an audit trail and is not
the accepted unified certificate.

`fullq3bcd0006_machine_c_prefix40.log` is the independent reverse-Pieri
diagnostic through degree 40.  Its 2,845 printed moments agree exactly with
the bounded-Littlewood implementation, and both implementations report 7,484
residual pairs with minimum `50`.  It is an implementation cross-check, not
the promotion certificate.  The final-source checker
`character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp` consumes
this prefix and the accepted bounded-Littlewood ledger.  Above degree 40 it
evaluates the determinants independently using ordinary modular arithmetic
and Newton forward differences at the integer nodes.  It checks exactly 538
higher moments modulo 21 primes whose 628-bit product exceeds the 607-bit
character bound, which makes every checked candidate unique, and then checks
the remaining 5,509 hierarchy values with exact GMP integers.  Its terminal
marker is `FULL_Q3_MODULAR_CHECKER VERIFICATION: ALL PASS`.
`fullq3bcd0007_machine_c_launch_failure.log`
records a wrapper failure before the verifier was invoked (`/usr/bin/time`
was absent); it is retained solely to make the run-ID history auditable.
