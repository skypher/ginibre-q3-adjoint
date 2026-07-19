# Full-Q3 finite classical residual run

## `fullq3bcd0001`

- Theorem slot: complete the exact residual boxes left after the certified
  tails for `B_18..B_21` and `D_31..D_70`.
- Host: machine C, `149.50.101.50` (`nb1cb2f`), 16 OpenMP workers.
- Arithmetic: C++20 with GMP `mpz_class`; no floating-point value enters a
  moment, hierarchy value, sign decision, or threshold comparison.
- Generator: `character_ring_iter/classical_boundary_certificate.cpp`, mode
  `--full-q3-b18-d31-residual-certificate`.
- Remote directory: `/root/q3_full_residual_20260715`.
- Remote log: `/root/q3_full_residual_20260715/fullq3bcd0001.log`.
- Command:

  ```text
  OMP_NUM_THREADS=16 ./classical_boundary_certificate \
    --full-q3-b18-d31-residual-certificate --progress
  ```

- Scope: all 137 higher-minus residual pairs in the four type-B rows and all
  4,732 pairs in the forty type-D rows.  The same run regenerates the stable
  moment prefix, every required finite-rank Pieri correction, and every
  hierarchy value.
- Kill condition: any moment-source mismatch, empty or inconsistent target
  ledger, nonpositive invariant-tensor moment, nonpositive hierarchy value,
  incorrect scope count, process failure, or resource exhaustion.
- Promotion target: a numbered exact-residual proposition in
  `full_q3_extension.tex`, shrinking the remaining classical obligation to
  `B_2..B_17`, `C_2..C_28`, and `D_4..D_30`.
- Merge/check: copy the complete log back under
  `certificates/full_q3/`, record its SHA-256, require the terminal marker
  `FULL_Q3_B18_D31_RESIDUAL VERIFICATION: ALL PASS`, and independently rerun
  the short directed-MPFR tail verifier before manuscript promotion.

### Result

Completed successfully on July 15, 2026.  The imported raw transcript is
`certificates/full_q3/fullq3bcd0001_machine_c.log`, with SHA-256
`82fd1d94a7799dcb57fbfb0e13677f7c00bd1449518cb92b72de019d10709c5e`.
It reports the exact scope `137 + 4732 = 4869`, no failed inequality, and
minimum

```text
D_31 a=8 n=15 value=1272988116891284109857142182380802.
```

The terminal success marker is present.  This first run is retained as an
audit trail; the source-matched promotion replay is `fullq3bcd0002` below.

## `fullq3bcd0002`

- Host: machine C, `149.50.101.50` (`nb1cb2f`), 16 OpenMP workers.
- Remote directory: `/root/q3_full_residual_20260715_v2`.
- Generator: `character_ring_iter/verify_full_q3_bd_residual_gmp.cpp`, with
  the included legacy engine and onset ledger authenticated independently.
- Recorded source SHA-256 values:

  ```text
  423d8aa6d23cff4ce8aaa765aa83421b6a6a63cc42aca45495c869a3fc88c38d  verify_full_q3_bd_residual_gmp.cpp
  bc22efcafafe2d40b5f5430cb8f9cb046bc6b9b188da684e61dc6fdfaa27077e  classical_boundary_certificate.cpp
  8d27088c9b9af35e7aca5eca5e99ddab70d38afdfb5e73b8fbd0c4cd1b7c6b04  full_q3_bcd_low_tail_data.hpp
  ```
- Imported transcript:
  `certificates/full_q3/fullq3bcd0002_machine_c.log`, SHA-256
  `9043a531c3411d36cb0aad1740274caaad5768ec7e059fbcb39593e7e78c8a14`.
- Acceptance: the transcript repeats the exact scope and minimum from
  `fullq3bcd0001`, ends with the required terminal success marker, and then
  records `__EXIT_STATUS=0`.

This source-matched replay meets the promotion target.  The result is
theorem-slot closure, not merely diagnostic evidence.

## `fullq3bcd0006_prefix40` (diagnostic only)

- Purpose: resource and arithmetic prefix audit for the later finite boxes
  `B_2..B_17`, `C_2..C_28`, and `D_4..D_30`; it is not a promotion replay.
- Host: machine C, `149.50.101.50` (`nb1cb2f`), AMD Ryzen 9 7950X,
  32 logical CPUs and 2 OpenMP workers.
- Start: `2026-07-15T04:29:59Z`.
- Remote directory: `/root/q3_full_remaining_20260715_v4`.
- Remote log: `fullq3bcd0006_prefix40.log`.
- Command scope: `--maximum 40 --progress`; consequently the program prints
  `full_certificate=0` and can end only with the diagnostic-prefix marker.
- Source SHA-256 values:

  ```text
  8de580680c75b05a21a9acba956babbf279813d6f553239ca4290028104f14c6  verify_full_q3_bcd_remaining_gmp.cpp
  bc22efcafafe2d40b5f5430cb8f9cb046bc6b9b188da684e61dc6fdfaa27077e  classical_boundary_certificate.cpp
  c3078f72b1cc02f22f7e6bb362ab9fad8e404b4854b319bc00d2775dd72592e6  full_q3_bcd_remaining_data.hpp
  bb7bcfd1869a16fe66fe68a0e15a9f9ae995916093e0f5aa562f69088f33e652  references/oeis_A002137_stable.txt
  ```

The run completed successfully at `2026-07-15T05:39:43Z`.  It checked 7,484
residual pairs through moment 40, with exact minimum `50` at `D_5,a=2,n=1`,
and ended with exit status zero.  These sources predate the final scope-count
and cap-domain assertions now present locally, so the run remains diagnostic.
Its 2,845 printed moment values were subsequently matched exactly by the
bounded-Littlewood implementation described below.

## `fullq3bcd0007_bounded_littlewood` (planned promotion replay)

- Theorem slot: close all 12,993 residual hierarchy inequalities in the
  remaining `B_2..B_17`, `C_2..C_28`, and `D_4..D_30` finite boxes.
- Structural mechanism: bounded Littlewood determinant identities reduce each
  finite-rank adjoint moment to a fixed-size truncated power-series determinant;
  evaluation over exact prime fields followed by GMP CRT reconstructs the
  integer moments.
- Host: machine C, `149.50.101.50` (`nb1cb2f`), 24 OpenMP workers.  The
  prime-field evaluations are independent parallel tasks; the CRT merge is in
  a fixed descending-prime order.
- Remote directory: `/root/q3_full_bounded_littlewood_20260715`.
- Remote transcript:
  `/root/q3_full_bounded_littlewood_20260715/fullq3bcd0007_machine_c.log`.
- Command:

  ```text
  OMP_NUM_THREADS=24 ./verify_full_q3_bcd_bounded_littlewood_gmp --progress
  ```

- Arithmetic: exact 32-bit Montgomery prime fields for determinant series and
  exact GMP integers for CRT reconstruction, character bounds, hierarchy
  values, sign tests, and scope counts.  No floating-point arithmetic enters a
  certificate claim.
- Source SHA-256 values before launch:

  ```text
  61bde65cd2ad346c1b128ac2d47b9353ccdf7cd1b043ada604eccf7cb380d656  verify_full_q3_bcd_bounded_littlewood_gmp.cpp
  7e913141d86e3f52d003b8cf5a1f2a0582d4a93722ecc60c04201b062acdabb9  full_q3_bcd_remaining_data.hpp
  ```

- Independent controls already completed: all 1,120 moments through degree 15
  and all 2,845 moments through degree 40 agree byte-for-byte (after label
  normalization) with the independent reverse-Pieri implementation; the latter
  comparison also agrees on all 7,484 hierarchy values' scope and minimum.
  A four-thread ASan/UBSan replay through degree 15 exits successfully.
- Kill condition: any source-hash mismatch, field exception, CRT bound failure,
  stable-prefix mismatch, `B_2/C_2` mismatch, nonpositive hierarchy value,
  scope other than 12,993, process failure, or missing terminal marker.
- Promotion target: replace the remaining-obligation remark by an exact
  bounded-Littlewood residual proposition and then state the all-classical and
  full adjoint-generated Q3 corollaries.
- Merge/check: import the immutable raw transcript under
  `certificates/full_q3/`, record source and transcript SHA-256 values, require
  `FULL_Q3_BCD_BOUNDED_LITTLEWOOD VERIFICATION: ALL PASS` and exit status zero,
  then run the separate exact/directed tail replay before promotion.

### Launch correction

Run `fullq3bcd0007_bounded_littlewood` exited with status 127 before invoking
the verifier because machine C does not install GNU `time` at `/usr/bin/time`.
Its short launch-failure transcript is retained and is not certificate
evidence.  The unchanged authenticated source and binary are relaunched as
`fullq3bcd0008_bounded_littlewood`, with the same host, 24-thread assignment,
kill condition, and promotion target; only the wrapper now uses Bash's `time`
keyword.  The promotion transcript is therefore
`fullq3bcd0008_machine_c.log`, never the `0007` launch record.

### Result of `fullq3bcd0008_bounded_littlewood`

The replay completed successfully on machine C from
`2026-07-15T05:47:19Z` to `2026-07-15T05:47:24Z` (4.423 seconds wall time,
68.539 seconds aggregate user time).  The imported immutable transcript is
`certificates/full_q3/fullq3bcd0008_machine_c.log`, with SHA-256

```text
0a5eadee6dd0fcde70c85774a34d274604cd2db9e17bff17eabf991282aadbb9.
```

It reports `maximum=59`, `parallel_threads=24`, reconstructs every requested
moment inside its exact character bound, checks the exact scope
`pairs=12993`, and finds minimum `50` at `D_5,a=2,n=1`.  The terminal success
marker and wrapper exit status zero are both present.  The source hashes in
the transcript match the pre-launch values above.  This meets the residual
promotion target; it is theorem-slot closure, not diagnostic evidence.

## `fullq3bcd0009_bounded_littlewood_and_exact_tails` (planned unified replay)

The `0008` transcript proves the complete residual box.  Before manuscript
promotion, the verifier was strengthened to consume the same reconstructed
moments in all six exact polynomial one-sided-tail ledgers and to rerun all
six exact rational Weyl-cap ledgers.  This makes the load-bearing exact part
of the remaining-row proof one source-matched replay rather than relying on
the old reverse-Pieri program's incomplete degree-40 tail prefix.

- Named theorem slot: the residual proposition together with the exact-tail
  proposition for `B_2..B_17`, `C_2..C_28`, and `D_4..D_30`.
- Host and command: machine C, 24 OpenMP workers,
  `OMP_NUM_THREADS=24 ./verify_full_q3_bcd_bounded_littlewood_gmp --progress`.
- Remote transcript: `/root/q3_full_bounded_littlewood_20260715/`\
  `fullq3bcd0009_machine_c.log`.
- Updated verifier SHA-256:

  ```text
  74632672afa53b8dbd99f7b1a0d6ebf10e7469dec6b41f0eae0cf7add5df1b23
  ```

- Additional kill conditions: any nonpositive exact polynomial input or
  margin, invalid Weyl-cap domain or geometry ledger, nonpositive cap margin,
  or exact-tail scope other than six polynomial plus six rational-cap rows.
- Promotion target: the same final all-classical and adjoint-generated Q3
  theorems, now with every exact tail dependency replayed in the accepted
  transcript.  The independently archived 58-row directed-MPFR transcript
  remains the complementary interval-tail certificate.

### Result of `fullq3bcd0009_bounded_littlewood_and_exact_tails`

The unified replay completed successfully on machine C from
`2026-07-15T05:57:25Z` to `2026-07-15T05:57:30Z` (4.471 seconds wall time,
70.360 seconds aggregate user time).  The imported transcript is
`certificates/full_q3/fullq3bcd0009_machine_c.log`, with SHA-256

```text
5efb579baad981d365deafa1ea15de3c8609db1b67b1465db5d2dee1b8c0b34c.
```

It authenticates the updated verifier and unchanged row ledger, asserts
`polynomial_rows=6` and `rational_cap_rows=6` after printing a positive exact
margin for every row, then repeats the complete residual result
`pairs=12993`, minimum `D_5,a=2,n=1,value=50`.  The terminal marker and exit
status zero are present.  This is the accepted exact promotion replay; the
earlier `0008` residual-only transcript is retained as supporting history.

## `fullq3bcd0010_fail_closed_final` (accepted final replay)

A strict post-replay code audit found no arithmetic discrepancy, but hardened
the verifier's fail-closed interface before final acceptance.  It now checks
the complete consecutive rank sequence, all 70 row coverage inequalities,
the exact method counts `6 + 6 + 58`, every rational-cap onset, and the maximum
moment endpoint before doing arithmetic.  It also prints each row's residual
count, allowing an independent sum back to 12,993.  OpenMP pragmas are guarded
so the same source compiles under both strict serial and strict OpenMP builds.

- Final verifier SHA-256:

  ```text
  257bcd0deaffea437ad25d2c76422b465ad33d7c2669df3d0d3b72c84134223d
  ```

- Pre-launch controls: strict serial and OpenMP builds pass; their 1,120
  degree-15 moments agree with each other and with reverse Pieri; a complete
  local replay reproduces every moment from `0009`, and the 70 printed row
  counts sum independently to 12,993.
- Host, assignment, command, kill conditions, and promotion target are
  unchanged from `0009`: machine C, 24 workers, full degree 59, all exact
  tails, residual scope, terminal marker, and exit status.  The promotion
  transcript is `fullq3bcd0010_machine_c.log`; `0009` is retained as the
  immediately preceding arithmetic-identical audit trail.

### Result of `fullq3bcd0010_fail_closed_final`

The final replay completed successfully on machine C from
`2026-07-15T06:15:10Z` to `2026-07-15T06:15:15Z` (5.221 seconds wall time,
85.999 seconds aggregate user time).  The imported accepted transcript is
`certificates/full_q3/fullq3bcd0010_machine_c.log`, SHA-256

```text
2f6d6673626efc2c834d36ecac5e02554aacc0ecf5a9478d3c1f6e7ae57d3291.
```

It authenticates the final source, checks all `6 + 6 + 58` method rows,
prints 70 consecutive row summaries whose pair counts independently sum to
12,993, reproduces every `0009` moment and exact-tail margin byte-for-byte,
reports minimum `D_5,a=2,n=1,value=50`, and ends with the required terminal
marker and exit status zero.  This is the final promotion transcript.

After importing the transcript, the final source was also compiled with
strict warnings plus AddressSanitizer and UndefinedBehaviorSanitizer and run
through the complete degree-59 calculation with 16 OpenMP workers.  It again
printed 70 rows summing to 12,993 pairs and the terminal pass marker, with no
sanitizer finding.  The local audit log SHA-256 was
`455c2f69497cfb1f5b9daa65d5757b825ae9d171a23560366a20806267613006`.

## Directed analytic replay and the type-B centering correction

The final source audit reread Courteaut--Johansson--Lambert's Proposition 3.3
at the random/full-trace boundary.  For `SO(2b+1)`, the random eigenvalue sum
has mean `-1` and the deterministic eigenvalue is `+1`; the full defining
character is therefore centered.  An interim wrapper run retained only the
deterministic contribution and is superseded.

- `fullq3bcdanalytic0001_machine_c_launch_failure.log` failed closed because
  the wrapper omitted `full_q3_bcd_low_tail_data.hpp`.  It has exit status one
  and is launch diagnostics only.
- `fullq3bcdanalytic0002_machine_c_superseded.log` compiled and ran the
  interim noncentered formula.  Its bytes remain authenticated, but its
  mathematical formula is wrong and it is not certificate evidence.
- `fullq3bcdanalytic0003_machine_c.log` is the accepted corrected replay on
  machine C (`nb1cb2f`, Ryzen 9 7950X, 32 logical CPUs).  It authenticates all
  four source/header inputs, records three successful compilations and three
  successful runs, checks all 768 finite middle-rank rows and all 44
  supplemental rows, and ends with both directed-MPFR terminal markers and
  aggregate exit status zero.  Its SHA-256 is
  `39f943dd18d4e6842d2dcd011ca0519edfc77f88579c4cc7ad58c854276b8ccf`.

The final middle verifier was additionally rebuilt with strict warnings,
AddressSanitizer, and UndefinedBehaviorSanitizer and completed all 812 rows
without a sanitizer finding.  The least supplemental `rho` is the centered
value `0.022694534162543144841871416180837` at `B_19`.
