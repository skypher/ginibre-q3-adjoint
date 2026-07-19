# Full-Q3 clean frontier fleet replay

## Runs `fullq3cleanfleet0002`, `fullq3fleetbmid0001`, and `fullq3fleetcheavy0004`

- Final-facing obligation: complete an independent recomputation of every long
  H25--H29 exact B/C frontier consumed by the companion Parts I--II proof, so
  that the final Part III import of the two-minus theorem has a current-source
  fleet audit rather than waiting for one serial host.
- Source relation: the final/current and older clean snapshots have 77 of 85
  replay-manifest inputs byte-identical.  All 27 long frontier C++ sources are
  identical.  The current `character_ring_iter/Makefile` changes only by adding
  the independent full-Q3 bounded-Littlewood target; existing frontier build
  rules and flags are unchanged.
- Computation: 71 independent exact C++/GMP tasks: 64 B shards and four C-only
  stages for H25--H28, followed by the H28 B absorption and the H29 C frontier
  and B absorption.  `run_full_q3_frontier_fleet_queue.py` checks the exact
  driver scope string, `SUMMARY failures=0`, process status, timeout, and source
  identity for every task.
- Final accepted partition: machine A (`optimus`, EPYC 7551) supplies the 24
  H25--H26 tasks from its four-worker/eight-thread queue; machine B supplies
  seven low-H27 tasks; machine C supplies the 40-task heavy remainder.  The
  exact disjoint selector is `run_full_q3_frontier_partition.py`; its three
  cardinalities are `24`, `7`, and `40`, and their union is checked against all
  71 task names produced by the original ledger.
- Run root on A: `/tmp/fullq3cleanfleet0001_a`.
- Orchestrator identity:
  `28b51756e053cb160e9cf5dde4f97a27ae416f084f7ef23ffedd02409098464a`.
  The staged copy is byte-identical to
  `run_full_q3_frontier_fleet_queue.py`; its 71 task names are unique and its
  seven C++ sources are byte-identical to the final workspace sources.
- Launch: `2026-07-15T09:22:15Z`, detached orchestrator PID `553826`, at
  niceness 10.  A one-stage C-only smoke replay passed before launch and was
  not included in the 71-task result directory.
- Because the original machine-A executor has all 71 tasks queued, the
  fail-closed watcher `monitor_full_q3_machine_a_low.py` freezes its process
  group immediately after all 24 H25--H26 PASS records are durable, checks the
  corresponding log hashes and both program/orchestrator zero-status footers,
  and writes only `machine_a_low_logs.sha256`.  Its SHA-256 is
  `59290ccecea9d51f401a7c7637587656bd80c69f57711a5f977107ecd0edd867`;
  partial later-task logs, if any, are excluded from the accepted manifest.
- Machine-B launch: `2026-07-15T10:10:30Z`, seven tasks, one process with 24
  OpenMP threads, run root `/root/fullq3_frontier_b_20260715`.
- At `2026-07-16T00:20:21Z`, during the `q=29` shard, machine B had about
  `3.6` GiB available and no swap while the authenticated trajectory still
  rose from `16780520` to `20720864` states.  To prevent an operating-system
  OOM kill, a dedicated 32-GiB `/swapfile_fullq3` was enabled before any page
  had been swapped.  This is a resource-safety measure only: it changes no
  source, binary, command, environment ledger, exact arithmetic, or accepted
  output.  Swap occupancy is monitored while the B queue remains live.
  At `2026-07-16T00:37:52Z` the shard had emitted the authenticated
  `after_removing_to=59` progress row with `20720864` states; available memory
  was about `825` MiB, the dedicated swap still contained zero bytes, and the
  kernel log contained no OOM record.  This is the peak-neighborhood snapshot,
  not an arithmetic PASS record; the task remains uncounted until its complete
  footer and partition digest are durable.
- Current machine-C candidate: launched `2026-07-15T23:05:56Z`, 40 tasks,
  one process with eight OpenMP threads, run root
  `/root/fullq3_frontier_c_20260715_w1`.
- Machine C is guarded independently by
  `monitor_full_q3_resource_envelope.py` (SHA-256
  `49ff73fe1eebbaf0fef5344f0d8884103fa28d57c1dfd00e6723f20a09143040`).
  It terminates and rejects the process group if available memory falls below
  12 GiB or swap use exceeds 2 GiB; its receipt is operational evidence only,
  never a substitute for any arithmetic task log.
- Kill condition: any source mismatch, missing binary, timeout, nonzero exit,
  missing exact scope header, missing zero-failure summary, duplicate task, or
  final coverage other than 71/71.  On failure, queued tasks are cancelled and
  no fleet summary is promoted.
- Merge/check: each accepted partition provides one log per task and a sorted
  SHA-256 manifest.  The final audit rejects duplicates, omissions, any union
  other than 71/71, any mismatch with clean stages 84--154, or any source,
  binary, scope, exit, and terminal-marker discrepancy.  It then combines the
  frontier ledger with the current changed-source stages and the H8--H24 clean
  transcript and compares the certificate union exactly with the paths
  reachable from the companion main theorem.
- A pre-promotion adversarial audit found that each fleet C++ program emits an
  internal `__EXIT_STATUS=0` and the Python orchestrator correctly appends a
  second process-status line.  The initial composition verifier incorrectly
  required only one occurrence and would have rejected valid logs.  The
  corrected verifier requires exactly two zero statuses, requires the second
  in the structured timing footer, binds every top-level PASS digest to the
  accepted task log, and independently compares the fleet command, shard,
  environment, and scope marker with the clean-driver stage.  Only the OpenMP
  thread count may differ across hosts.  Replay PASS, exit, build, wrapper, and
  partition-terminal records are unique exact lines, not substrings.  It
  additionally requires the complete
  ordered arithmetic projection of each task to equal its uniquely scoped
  block in the accepted paper-cited H25--H29 transcript; omission, insertion,
  reordering, and alteration negative tests all fail.  Exact rank-indexed
  output cardinalities are checked independently.  It also requires the seven
  cited frontier certificates to contain exactly the 71 ledger scopes, so an
  extra stale or unreplayed certificate block is fatal.  The same composition now
  binds the 34 direct certificate replays and all 19 H23--H24 scoped blocks,
  for 124 exact certificate-data projections in total; declaring a certificate
  name without matching its data is insufficient.
- Promotion target: unconditional current-source verification of the companion
  two-minus theorem dependency used by `full_q3_extension.tex`.

The promotion collector binds machine A's copied top transcript to the exact
SHA-256 stored by the freeze watcher, rather than trusting the terminal marker
alone.  It also parses every archived machine-C resource-watch line, rejects
malformed records, and independently enforces the 12-GiB available-memory
floor and 2-GiB swap-use ceiling before admitting the arithmetic evidence.

### Rejected resource-envelope launches

The initial eight-worker launch began at `2026-07-15T08:40:50Z`.  Before any
task completed, the eight H25 processes reached approximately 90 GiB aggregate
RSS and the highest shards were still growing.  At `2026-07-15T09:20Z` the
orchestrator PID `528075` and its eight explicit children were terminated as a
resource-safety failure.  Its eight incomplete logs and top-level transcript
are preserved under `results_aborted_workers8_20260715T0920Z` and
`fullq3cleanfleet0001_a_aborted_workers8.log`; none is admissible evidence.
Run `fullq3cleanfleet0002` changes no arithmetic source or task boundary.  It
adds an explicit run identifier and uses four workers to keep the worst
simultaneous GMP state maps within machine A's physical-memory envelope.

The first machine-C heavy-partition candidate `fullq3fleetcheavy0001` used
four workers.  Before any task completed, its aggregate task RSS grew from
about 21 GiB to 49 GiB in four minutes and was still increasing.  The four
explicit children and their orchestrator were terminated at
`2026-07-15T10:33Z`; partial logs are preserved under
`results_aborted_workers4_20260715T1033Z` and are not evidence.  Candidate
`fullq3fleetcheavy0002` used three workers.  Its three first H27 ranks grew to
approximately 95 GiB aggregate RSS, left less than 1 GiB available, and forced
7.2 GiB into swap before any task completed.  It was terminated at
`2026-07-15T11:13Z`; its partial logs are preserved under
`results_aborted_workers3_20260715T1112Z` and are not evidence.  Candidate
`fullq3fleetcheavy0003` used two workers and completed eight tasks before the
independent monitor measured `12381884` KiB available, below its fail-closed
`12582912`-KiB floor, at `2026-07-15T16:46:10Z`.  The monitor terminated the
process group and wrote an aborted receipt; the complete candidate directory
is preserved at `/root/fullq3_frontier_c_20260715` and none of its logs is
accepted into the final partition.  Current run `fullq3fleetcheavy0004`
changes no source, task selection, or arithmetic and uses one worker; it
becomes accepted evidence only after all 40 tasks and the final composite
verifier pass.
