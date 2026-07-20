# Clean-room replay contract

## Command

The full contract is an expensive exact replay.  Run it from the repository
root on machine C (the default primary compute host):

```text
make -C ginibre_q3 clean-room-replay
```

CPU affinity, available host RAM, and cgroup memory limits are detected by
default.  An explicit value remains a hard upper bound:

```text
make -C ginibre_q3 clean-room-replay REPLAY_THREADS=32
```

Success is the final line

```text
[replay] ALL GINIBRE Q3 CHECKS PASSED
```

and a zero process exit status.  Any missing input, hash mismatch, nonzero
accepted replay, build failure, failed inequality, timeout, or unresolved
LaTeX reference makes the command fail.

Before an expensive replay, run the complete publication preflight:

```text
make -C ginibre_q3 publication-preflight
```

It runs the document and dependency checks, moment and cone formulas,
classification coverage, the exact `SO(3)` obstruction, the cross-document
two-minus interface, the independent frontier-formula audit, and every
artifact manifest.  These are structural, formula, and byte-authentication
checks; they do not re-evaluate the load-bearing arithmetic signs.

## Verification levels and trusted base

The archive deliberately separates three claims that had previously been
described together.

1. `publication-preflight` checks theorem/proof structure, the explicit
   residual B/C caller contract, classification coverage, independent small
   formulas, and manifest integrity.  It is suitable for an ordinary
   workstation and completes without regenerating large moment tables.
2. An accepted transcript authenticates one completed exact or
   directed-rounding execution.  Its hash proves byte identity only; reading
   a positive margin from that transcript is not an independent
   recomputation of the margin.
3. `clean-room-replay` rebuilds the active C++ programs and re-evaluates every
   sign used by Parts I--II.  This is the arithmetic proof replay and is the
   only supplied command that independently recomputes all such signs.

PDF synchronization is a separate fail-closed level.  The command
`make -C ginibre_q3 publication-artifact-audit` performs three LaTeX passes
on each formal component, rebuilds the reader PDF from fresh Parts I--II and
compact Part III components, and verifies all four resulting SHA-256 values against
`publication_artifacts.sha256`.  It detects a checked-in PDF made from older
source even when that stale PDF is nonempty and valid LaTeX output.

The trusted software base for level 3 is the manifested C++ source, the C++
compiler and standard library, GMP, MPFR, OpenMP, `pdflatex`, and the Python
orchestration that checks stage status and coverage.  The submission does not
claim formal verification of this software.  In particular, the short suite
must never be described as a lightweight checker for all arithmetic
certificates; it checks their scope and authentication, not their numerical
contents.

`ENVIRONMENT.md` records one fully specified validated software stack and the
host on which the July 2026 publication checks were rebuilt. The Makefiles
remain authoritative for flags, while that environment record removes any
ambiguity about the tested compiler, GMP/MPFR, OpenMP, Python, and TeX
versions.

The exceptional-prefix checker is fail-closed: it rejects an input that ends
before `m_5`, a gapped moment sequence, conflicting duplicate moment rows, a
negative `Q_3`, or a negative Chain difference, and every such rejection has
a nonzero process status.  The clean-room driver separately requires the
literal terminal marker `Chain verification: ALL PASS`.  The failure paths
are regression-checked with short synthetic moment streams; all five
authenticated G2/F4/E6/E7/E8 prefixes retain their byte-identical successful
output after this hardening.

The active type-B H8--H27 frontier is a single 337-case hybrid exact
certificate: hook-length sums close 296 cases analytically, and a ten-prime
bounded-Littlewood determinant/CRT reconstruction closes the remaining 41
low-rank cases through moment 67.  It uses the stable moment as a proved
uniqueness bound and replaces all mandatory reverse-Pieri frontier shards.
The 337 independent paired-row partition sums and the prime-by-interpolation-
node determinant grid both use the affinity/RAM-limited OpenMP team; the
ordered reduction remains deterministic.  This stage takes 24.93 seconds on
four cores, down from 35.45 seconds after parallelizing only the grid.
The H9--H27 mixed programs now run only their
closed type-C formulas; archived B shards remain optional controls.  The
all-range half-stable bridge and power-loss verifier remain exact GMP scans.
The clean build likewise compiles only executables invoked by the mandatory
replay; seven optional-control binaries remain available as explicit Make
targets but are no longer paid for on the four-minute critical path.  The 19
closed-form frontier executables share one unoptimized unity translation unit
and dispatch by executable name; they still run in hundredths of a second.
Optimization remains enabled for every determinant, recurrence, MPFR, and
character-ring kernel.
The half-stable bridge is capped at 32 OpenMP workers because it keeps a large
exact state map per worker.  Every executable clean-room stage now has a
hard four-minute ceiling (`REPLAY_MAX_STAGE_SECONDS=240`).  A stage that
crosses it is terminated with its isolated workspace and log preserved; it
must be profiled and optimized or replaced by a proved analytic reduction
before it may re-enter the replay.

After the one-time build and stable-recurrence check, the compute-dense active
type-B frontier receives the full adaptive core allocation once.  The four
remaining dependency-disjoint proof groups then run concurrently.  The
resource planner caps compiler jobs at
eight, reduces simultaneous proof groups only below the audited RAM floor,
and assigns explicit endpoint, B/C, exceptional, and classical CPU shares.
After the exact bridge, Weyl-orbit, and adaptive-grid reductions, a four-core
host exposes four endpoint workers, two exceptional workers, and one worker
in each short B/C and classical group.  The brief initial overlap is
time-shared; as the two short groups drain, their cores are automatically
reclaimed by the two long-lived groups instead of sitting idle.  On the
validated 64-core host it selects 8 endpoint
threads, 24 B/C-residual threads, 24 exceptional threads, and 8 classical-tail
threads.  Measured kernel saturation caps prevent a larger host from adding
counterproductive workers.
Per-stage logs, required-marker checks, certificate accounting, and the
240-second ceiling remain independent and fail closed under this scheduling.
The driver also enforces a 240-second budget for the entire isolated workflow,
including integrity audits, compilation, arithmetic, and the compact PDF
rebuild (`REPLAY_MAX_TOTAL_SECONDS=240`).
The 64-core, 112.7-GiB clean-room run of 2026-07-20 completed all 105 stages,
71 manifests, 3,334 hashes, and the 53-page PDF in 168.97 seconds.  The
fixed-frontier bridge itself took 38.36 seconds under full replay load, down
from 88.69 seconds before this reduction.
The conservative fixed-boundary/log-convex bridge passed a four-core
standalone replay in 17.45 seconds.  Extending its shared audit through index
30,898 with rigorous 192-bit dyadic prefixes gives a complete one-core bridge
time of 21.90 seconds.  Adaptive exact E8 meshes reduce that
four-core stage from 83.14 to 22.23 seconds, while signed-permutation orbit
reduction cuts the complete low-rank Weyl source from 68.70 to 2.37 seconds.
The identity
`binom(j,q)*s[j-q]/s[j] = (s[j-q]/(j-q)!)/(q!*s[j]/j!)`, together with the
same exact log-convexity audit through index 30,898, reduces the power-loss
stage from 3,904,626 comparisons to 402 boundary comparisons; it takes 1.27
seconds on one core.
The 21 finite type-D source rows are now reconstructed together by one exact
bounded-Littlewood determinant/CRT grid through the largest consumed degree
46.  It checks 1,060 authenticated moment/correction claims in about five
seconds on four cores, replacing the repeated rank-by-rank Racah--Speiser
growth; the separate D4 modular/Weyl control is retained.

The exceptional source audit applies the same ceiling to character-ring
growth.  Fresh Racah--Speiser regeneration stops at `G2:m38`, `F4:m65`, and
`E6/E7/E8:m42`; these prefixes validate the Cartan data, Weyl reflections,
and character-pairing implementation.  Every longer moment actually consumed
by the proof is then compared exactly with the corresponding `g[2]`, `f[4]`,
`e[6]`, `e[7]`, or `e[8]` block in the BPV ancillary source, after which the
full GMP Chain checks run on the complete consumed ledgers.  The old attempt
to regenerate the archival `F4:m220`, `E6:m80`, `E7:m70`, and `E8:m100`
suffixes in one process is not a replay obligation: its support grows to tens
of millions of highest weights and individual historical steps exceed five
minutes.

Success verifies the certificate and artifact contract for the paper's
unconditional, exactly-two-minus repeated-adjoint theorem.  This is not a
claim about Ginibre's full unequal-character Q3 cone.  The old type-D envelope
`DEnv` is false and is not an allowed
condition.  The replay verifies its exact counterexample and the complete
unconditional replacement for `D_4..D_295`; no reachable conditional result
is permitted.

## Part III full adjoint-generated replay

The separate Part III theorem concerns every sign pattern on the closed
multiplicative cone generated by the adjoint character.  Its aggregate audit
and document build are:

```text
make -C ginibre_q3 full-q3-extension
```

The final finite classical stage is
`full-q3-bcd-bounded-littlewood-audit`.  It reconstructs all required
`B_2..B_21`, `C_2..C_28`, and `D_4..D_70` moments by exact bounded-Littlewood
determinants over prime fields, performs a GMP CRT merge beyond the character
bound, reruns the six exact polynomial and six exact rational-cap tails, and
checks exactly 17,862 residual hierarchy values.  This includes the formerly
separate 4,869-case B/D reverse-Pieri box.  The separate 58-row
directed-MPFR tail stage remains part of the aggregate target.  The accepted
unified current-source exact transcript is
`certificates/full_q3/fullq3bcdboundedfinal0001_current_source.log`; success
requires its verifier/data hashes, terminal marker, wrapper exit status zero,
and adjacent SHA-256 sidecar.  The older `fullq3bcd0010` and
`fullq3bcd0002` transcripts retain independent historical coverage of the
12,993-case and 4,869-case subledgers.  The mandatory modular/Newton result is
separately authenticated by
`certificates/full_q3/fullq3bcdmodularfinal0001_current_source.log` and its
adjacent sidecar.

The finite-middle and high-rank analytic branches are authenticated by
`certificates/full_q3/fullq3bcdanalytic0003_machine_c.log`.  In the type-B
Fredholm formula the fixed `+1` eigenvalue cancels the `-1` mean of the random
eigenvalue sum, so the full defining trace uses the centered Gaussian factor.
The accepted transcript authenticates that corrected source and the low-tail
header, checks 768 middle rows and 44 supplemental rows with 384-bit outward
rounding, and records all compile/run statuses as zero.  The earlier `0002`
analytic transcript is explicitly superseded; `0001` is a fail-closed launch
record.

Before arithmetic, `full-q3-document-audit` requires 51 numbered results and
51 adjacent proofs, unique labels, no missing references, exact citation-key
coverage, the SO(3) obstruction, and the centered type-B cancellation.
`full-q3-artifact-audit` verifies every repository manifest, including the
Part III source manifest.  A mismatch prevents the aggregate target from
building the PDF.

This proves the full Q3 quantifiers only on the adjoint-generated cone.  The
literal nonabelian analogue using every real positive-definite function is
false, with exact `SO(3)` value `-8/279` in Part III.

## Reproduction tiers and resource envelope

The archive separates checking from regeneration so that a reader can see
which conclusion each command supports.

| Tier | Command or target | Purpose | Documented resource scale |
|---|---|---|---|
| Preflight | `make -C ginibre_q3 publication-preflight` | Parse theorem interfaces and the explicit B/C contract, check formulas and coverage, and authenticate manifests; does not recompute all arithmetic signs | Seconds to a few minutes; ordinary workstation |
| Part III arithmetic | `make -C ginibre_q3 full-q3-extension` | Rebuild every Part III exact/MPFR verifier and its 35-page formal/computational supplement | Parallel exact suppliers; each underlying computation is required to remain below five minutes on the validated 64-thread host |
| Parts I--II regeneration | `make -C ginibre_q3 clean-room-replay` | Rebuild every main-theorem-reachable two-minus certificate in an isolated tree | Hard 300-second ceiling for each executable stage; a ceiling failure preserves its profiling log |

The Part III target verifies the imported two-minus theorem's source and
interface bindings but does not rerun its 200-stage arithmetic.  A complete
reproduction therefore requires both aggregate commands.  The files under
`referee_audits/` are author self-audits and redundancy checks, not independent
peer review.  SHA-256 records provenance and byte identity; only the exact or
outward-rounded comparisons specified in the numbered results decide signs.

## Final-source archival replay

The publication archive is produced only after a clean-commit replay of both
arithmetic aggregates. Choose an absolute log path outside the repository:

```text
FINAL_REPLAY_LOG=/absolute/path/final-replay.log \
  make -C ginibre_q3 final-publication-replay
```

`run_final_publication_replay.sh` rejects a dirty worktree and an existing log
path. It records the exact commit, final-source-manifest hash, deterministic
PDF epoch, host, compiler, Python version, and UTC start/pass time of every
stage. It uses unbuffered Python and the progress-enabled C++ entry points.
The terminal marker is `FINAL_PUBLICATION_REPLAY: ALL PASS`. It then verifies
that the replay did not alter the source tree.

After that marker is present, build the immutable-deposit payload with:

```text
make -C ginibre_q3 final-release-archive \
  FINAL_REPLAY_LOG=/absolute/path/final-replay.log \
  FINAL_RELEASE_ARCHIVE=/absolute/path/ginibre-q3-final.tar.gz
```

The archive builder again requires a clean worktree, refuses to overwrite an
output, includes every tracked `ginibre_q3` file plus the replay log and exact
commit metadata, normalizes ownership and timestamps, and prints the archive
SHA-256. Running it twice from the same commit and replay log produces the
same bytes. The resulting archive and checksum are the objects to deposit in
the public DOI-bearing repository; local tooling cannot assign that external
identifier.

## Compressed full-range implementation cross-check

The accepted bounded-Littlewood calculation and the reverse-Pieri traversal
are structurally independent algorithms. The archived reverse-Pieri prefix
checks all 2,845 required moments through degree 40, covering 7,484 residual
pairs. Extending that traversal through degree 59 creates a much larger GMP
frontier but adds only 538 moments needed by the remaining 5,509 pairs. The
fail-closed target therefore verifies precisely that smaller complement:

```text
make -C ginibre_q3 full-q3-bcd-independent-audit
```

The target builds
`character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp`. It first
requires exact equality with the archived reverse-Pieri prefix. For degrees
41--59 it evaluates the bounded-Littlewood determinants in ordinary modular
arithmetic, using Newton forward differences instead of the promotion
verifier's Montgomery arithmetic, Lagrange weights, and GMP CRT
reconstruction. All 538 candidate moments must agree modulo 21 distinct
primes. Their 628-bit product is larger than the independently enforced
607-bit maximum character bound; since both the true moment and candidate lie
between zero and that bound, the modular agreement proves equality over the
integers. The checker then evaluates all 5,509 remaining hierarchy values
with exact GMP integers and requires strict positivity, in addition to the
7,484 prefix count. It emits
`FULL_Q3_MODULAR_CHECKER VERIFICATION: ALL PASS` only after all scopes and the
modulus-dominance test pass.

The modular checker is a mandatory dependency of
`make -C ginibre_q3 full-q3-extension`; a successful Part III aggregate must
therefore pass both the promotion verifier and this distinct implementation.
The optional aggregate
`make -C ginibre_q3 full-q3-extension-independent-controls` additionally
runs the historical reverse-Pieri supplier for the 4,869-case B/D box.  That
larger traversal is retained as a further algorithmically independent overlap,
not as a second mandatory supplier.

To preserve a timestamped raw audit log, use the wrapper with a new absolute
directory:

```text
FULL_Q3_INDEPENDENT_LOG_DIR=/absolute/new/directory \
  ./ginibre_q3/run_full_q3_bcd_independent_audit.sh
```

The compressed target is included in `full-q3-extension`; it remains a
separate implementation and source file from the determinant promotion
verifier.

## Isolation guarantee

The driver creates a new temporary directory and copies only the Ginibre Q3
subproject into it.  Before compiling, it:

1. rejects symlinks and any manifest target resolving outside the copied
   tree, then verifies every `.sha256` manifest there (including duplicate
   entry rejection within each manifest);
2. checks every replay marked `accepted` in the post-`m=29` classification
   has `__EXIT_STATUS=0`;
3. removes every copied active ELF executable, Python bytecode cache, and
   main/full LaTeX output before rebuilding.  Checksum-bound ELF files inside
   inert historical certificate snapshots are preserved for manifest
   verification but are never invoked.

It then builds the active replayers from source.  Thus a tracked or untracked
binary in the working tree cannot satisfy an executable replay stage.  A
successful temporary tree is deleted; a failed tree and its per-stage logs
are preserved and printed for diagnosis.

## Recomputed proof spine

The command recomputes:

- the formal `SU(3)`, `SU(4)`, and `SU(5)` moment recurrences by exact
  Bessel-determinant expansion, binomial conjugation in the Weyl algebra,
  and differential reduction modulo `t*y''+y'-y=0`; the checker uses only
  integer Laurent-polynomial arithmetic and requires a zero residual;
- the SU(N) finite transition strip for `N=6..18`, the uniform strip base
  margin and symbolic monotonicity polynomials, and the endpoint/overlap
  certificates using exact GMP arithmetic;
- every type-D adjoint-moment source and bounded exceptional prefix routed to
  the active root-datum checker, by reconstructing its root datum, iterating
  the exact Racah--Speiser character ring, and applying character
  orthogonality; the 182 consumed low-rank B/C source claims are instead
  compared termwise with the bounded-Littlewood determinant reconstruction,
  while every longer consumed exceptional row is compared termwise with the
  BPV ancillary table before the full integer Chain audits run;
- the stable classical sequence from the proved exact three-term recurrence;
  the archived OEIS row is checked in every degree before any downstream
  correction stage, so each later table lookup is certified equal to the
  recurrence value;
- 57 classical `m=29` source rows by compiling the C++/GMP Weyl/Pieri
  generator: `B_2,C_2,B_3,C_3`, `B_19..B_30`, `C_19..C_30`, and
  `D_35..D_63`;
- 832 finite `B/C` corrections from a common modular Pieri traversal with
  exact GMP CRT reconstruction, including an exact comparison with all 182
  consumed archived claims in 14 low-rank rows, followed by all 870 row-gated
  Chain inequalities (`416` exact and `454` interval-minimized);
- the 57 first-hit trace rows, the 118 high-edge Chernoff rows, and the
  Rains-square cutoff for every integer `296 <= C <= 10000`, all in the
  same 384-bit directed-MPFR C++ verifier;
- the exact C++/GMP `D_4,j=18` counterexample to `DEnv`;
- the D4 modular reconstruction audit, the identity-checked machine-C exact
  adjoint-moment prefixes, and all 530 row-gated `D_4..D_24` Chain steps in
  exact C++/GMP interval arithmetic;
- all 36 nonvacuous row-gated steps for `D_53..D_63` from the proved
  nonnegative correction box, again in exact C++/GMP interval arithmetic;
- the exact-rational C++/GMP `D_10` cap/Chebyshev tail from odd exponent 63,
  the directed-MPFR exact-moment `D_11` MGF tail from odd exponent 75 and its
  exact GMP bridge through that onset, and the 13-row directed-MPFR
  exact-moment MGF tail for `D_12..D_24` from exponent 63;
- the type-B unitary-square, B/C tilted-MGF, and combined 333-row layered-MGF
  certificate, including 243 unconditional type-D onsets `D_53..D_295`;
- the 28 exact/mixed finite-rank type-D MGF onsets `D_25..D_52`, with exact
  SO-even defining MGFs for `D_25..D_28` and exact Rains square MGFs for all
  28 rows;
- all 707 correction-aware finite Chain steps for `D_25..D_52` using exact
  GMP two-sided correction boxes and OpenMP;
- an independent `D_19` exact finite-rank push-moment, directed-MPFR onset,
  and 24-step correction-aware GMP replay;
- the exact GMP A-free Chain frontier from those onsets across every
  `D_53..D_295` row and `m=30..163`;
- the corrected root-normalized caps for `B_4,B_5,B_6,B_8`, the 22-row
  C++/GMP/MPFR exact determinant-MGF tail for
  `B_7,B_9..B_13,C_4..C_19`, and the short exact-integer 34-step bridge for
  its six delayed onsets; the normalization-defective low-rank tail logs are
  excluded from the accepted boundary;
- the active 530-row B/C first-hit replay; the former type-D portion of
  the historical 796-row replay is replaced by the exact D25--D52 and
  A-free D53--D295 C++/GMP bridge verifiers above;
- the active quantified caller contract for the 58 family-wide residual B/C
  correction propositions at offsets `0..27`; this domain-aware audit checks
  their rank ranges, the exact table labels, tilted-tail onsets, and all-row
  GMP scope.  It independently proves that main-theorem propagation consumes
  at most offset `27` and enumerates the 2,834 high-rank odd targets that
  would be uncovered if the old rank-61 cap reappeared.  The additional
  offset-28 and offset-29 results are historical overlaps, not premises of the closure;
- the degree-8 cutoff-80 B/C Chebyshev negative-tail bound, the 402-row
  directed-MPFR interval onsets, the exact half-stable bridge through
  `m=15447`, and all `3,904,626` exact power-loss inequalities;
- the hybrid exact 337-case H8--H27 type-B hook-length/determinant certificate,
  the `B_14,j=37` seed, and the inexpensive type-C formula portions of the
  ninth through twenty-seventh frontiers; H28--H29 and the old type-B
  reverse-Pieri shards are authenticated historical overlaps;
- the G2/F4/E6/E7/E8 exact prefix audits after their source moments have been
  regenerated independently;
- the F4/E6/E7/E8 Dunkl coefficient identities using parallel GMP;
- the G2/F4 exact GMP rectangular tails, the rank-corrected E6/E7 OpenMP/GMP
  rectangular tails (including the E7 degree-26 exact moment majorant), every
  E8 negative-bound row through the parallel C++/GMP
  source/majorant/moment audit,
  every E8 rectangular bridge/tail row through exact GMP geometry with
  384-bit outward-rounded MPFR accumulation, and
  the E8 `m_0..m_100` finite bridge;
- three clean `pdflatex` passes for the self-contained Parts I--II manuscript
  `paper.tex`, and the absence of undefined references or citations;
- the formal document contract: `paper.pdf` is nonempty and self-contained;
  `full_q3_main.pdf` states the compact Part III proof spine; and the separately
  audited `full_q3_extension.pdf` is its load-bearing supplement.  The first
  two are collected in the 59-page `submission.pdf`; `paper_full.pdf` is an
  optional derivation archive rather than a formal submission component.

The clean-room arithmetic replay builds `paper.pdf` in its isolated tree but
does not publish that temporary file.  The final release driver therefore
runs `publication-artifact-audit` in the final-source tree after both
arithmetic aggregates.  This ordering prevents the temporary Parts I--II PDF
from being discarded while an older tracked `paper.pdf` is silently used by
`submission.tex`.

The replay also validates all archived SHA-256 manifests and
the accepted/diagnostic classification boundary before consuming any
certificate.  It enumerates every in-tree file argument consumed by the proof
spine and rejects any one absent from the verified manifests.  It additionally
parses the paper's identity ledger and rejects any certified-output or inline
SHA-256 that is not present in a verified manifest.  The reference audit
separately requires exact manifest coverage of
the whole `references/` directory and an `active_paper_sources.tsv` mapping
from every bibliography key to an archived file and verified source locus.
The active Etingof and Garibaldi--Guralnick--Rains mappings point to their
vendored primary PDFs, not to a secondary Markdown audit.
The source-structure audit pairs all numbered theorem-like environments with
their proofs, handles the intentionally deferred main-theorem proof, walks
the explicit cross-reference graph, and rejects enumerated unresolved-status
phrases reachable from the final assembly.  This is a syntactic consistency
check, not a theorem prover; mathematical dependency discharge is the content
of the numbered proofs in Parts I and II.
After all executable stages finish, the driver extracts every certificate log
cited by a main-theorem-reachable statement or proof and requires exact
equality with the set actually recomputed.  Every reachable log must be
hash-authenticated, and every reachable post-`m=29` log must also be classified
`accepted`.  The current proof spine contains 43 logs: 42 post-`m=29` logs and
the classical-trace first-hit log.  Adding a new reachable log without a replay
stage, retaining a replay after its proof dependency is removed, or exposing a
conditional or diagnostic post-`m=29` certificate to the main theorem makes
the command fail.

## Source-regeneration boundary

The archive contains raw character-ring, determinant-isotypic-sum, partitioned
frontier, and exploratory search transcripts.  A SHA-256 digest authenticates
only the bytes.  The standard command independently regenerates every
low classical character-ring claim and each bounded exceptional prefix passed
to a main-theorem source-replay stage.  Longer exceptional suffixes are
checked termwise against the independently archived BPV ancillary table; the
full downstream Chain calculation still consumes and checks every required
integer row.  The command rebuilds the theorem-reachable determinant/frontier certificates, regenerates
the stable classical moments from the proved recurrence while checking the
OEIS comparison table, verifies recorded exit statuses and the accepted
classification, and checks every downstream exact inequality.

The finite row-gated bridge through `m=29` is source-generated.  The command
reruns the common Pieri/CRT verifier for every `B/C` row and the exact type-D
interval verifier for the remaining high-rank rows; the active `D_4..D_52`
replays supply the other type-D steps.  The former 297-file bridge archive,
its classifications, and its transcript remain byte-authenticated provenance,
but no theorem consumes their statuses or printed margins.

In historical type-D lower-bound logs, a bare `D_b Delta_j = ...` line above
`j=2b` is only the determinant component, not the full correction.  It is
ignored as an exact moment by both the source-regeneration audit and the
downstream interval checker.  In-window bare lines and explicit full or
exact-Weyl moment lines are compared exactly.

Thus the one-command contract replays every main-theorem-reachable certificate
verifier, including the entire finite row-gated bridge, all 42 reachable
post-`m=29` certificate logs, and the reachable classical-trace first-hit log.
Historical exploratory searches remain outside the proof spine.

## Arithmetic policy

All load-bearing discrete arithmetic uses GMP integers or rationals.  Claims
involving logarithms, gamma functions, or analytic interval comparisons use
outward-rounded MPFR at the precision stated by the program.  OpenMP is used
for the partition, Dunkl, type-D frontier, interval, and E8 workloads.  Python
only orchestrates or performs short replay/audit checks in this command; no
expensive search or load-bearing certificate generation runs in Python.
Some exact-GMP programs additionally print ordinary floating-point `log10`
summaries to identify a worst row.  Those values are diagnostics, may differ
in their last displayed digits across hardware, and never select a certificate
sign or discharge an inequality; the corresponding pass/fail decision is the
exact integer comparison printed immediately before the summary.
