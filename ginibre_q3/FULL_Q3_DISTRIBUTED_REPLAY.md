# Distributed clean replay for the full Ginibre-Q3 proof

## Run `fullq3distributed0001`

This run composes the companion Parts I--II replay without weakening the
stage ledger of `clean_room_replay.py`.  The current driver expands to 200
fresh-build and arithmetic stages before its six final LaTeX passes.  The
composition is accepted only if `verify_full_q3_distributed_replay.py`
reconstructs those stages from the driver and finds exact, duplicate-free
coverage as follows.

| Clean stages | Count | Evidence |
|---|---:|---|
| 1--49 | 49 | final-source machine-C clean build and prefix through the B14 bad-shape stage |
| 50--83 | 34 | older machine-C clean snapshot, H8--H24, after byte-equivalence and build-rule audit |
| 84--107 | 24 | current-source machine-A partition, all H25--H26 tasks |
| 108--113, 129 | 7 | current-source machine-B partition, H27 B-shards through `q=32` and the H27 C task |
| 114--128, 130--154 | 40 | current-source machine-C heavy partition, the remaining H27 tasks and all H28--H29 tasks |
| 155--184 | 30 | current-source classical post-frontier replay on `optimus` |
| 185--200 | 16 | current-source exceptional replay on `optimus` |

The final-source machine-C prefix is frozen under
`certificates/full_q3/distributed0001/current_prefix`.  Its top transcript has
SHA-256 `bf9639bd614a12147af7ee1ae7a0604717707feffb38b537baa26d0a723da9b6`.
The load-bearing 15,417-case bridge log has SHA-256
`da663c330570cf1efbb7e30d4bc4d096a74b1ffbd6e12e0ba577d1f40178830d`,
and stage 49 has SHA-256
`a17a760d7219bdf0d4c5b4fc010f5ea92a5c55f2c55d5b69538a7cbf91e6b7e3`.
An independent invocation of the composite verifier's prefix routines found
49 distinct logs, all nine preflight markers, exact command/environment
agreement, and no failures.

The verifier checks the exact command, environment override, required scope
header, zero-failure marker, program-success and orchestrator-success evidence,
source identity, build identity, and log identity for every stage.  For fleet
stages it independently compares the pinned task ledger with the clean-driver
stage specification, permitting only the host-specific OpenMP thread count,
and requires the SHA-256 printed in each top-level PASS record to equal the
accepted task log.  The three frontier selections
have cardinalities `24+7+40=71`; their task names are unique and their union is
required to equal the original unpartitioned 71-task ledger.  It then requires
exact equality between the union of replay-declared certificates and all 43
certificate logs reachable from the companion main theorem.
All replay PASS, exit-status, wrapper-identity, partition-terminal, and build
identity records are required as unique complete lines; embedded substrings,
suffixes, malformed timings, and duplicate metadata are rejected.

Certificate-name coverage is not treated as certificate-data verification.
The verifier additionally compares 124 complete regenerated payloads with the
paper-cited data: 34 direct certificate replays, all 19 uniquely scoped
H23--H24 B/C blocks, and all 71 H25--H29 fleet blocks.  Equality is ordered and
complete.  The only normalizations are absolute path relocation, host thread
counts, and redundant scheduling-progress lines whose row data are retained
in the deterministic ordered tables.  Missing, inserted, reordered, or
numerically altered load-bearing lines fail the composition.  For the frontier
certificates, the verifier also requires equality between the complete set of
certificate scopes and the 71-task ledger: all seven cited certificate files
must contain exactly 71 unique blocks, with no stale or unreplayed extra block.

For H25--H29 the audit splits the authenticated certificate transcript cited
in `paper.tex` into its uniquely scoped task blocks, normalizing only the
host-dependent terminal OpenMP count.  After removing orchestration,
scope-header, and summary lines, the complete ordered list of regenerated
arithmetic lines must equal the corresponding accepted block byte for byte.
It separately requires the exact B-rank set in every shard, all B
count/stable/margin records, the complete C fixed-point-free and arithmetic
rank sets, and the eight- and ten-rank H28/H29 absorption ledgers.  Thus a zero
failure counter cannot mask an omitted, added, reordered, or altered datum.

The verifier also compiles the actual H25 reverse-Pieri predecessor routine
against an independent horizontal-two-strip interlacing oracle.  It checks all
1,595 partitions through size 18, verifies that the H25--H28 kernel source
blocks are identical, and checks the live encoding bound `121 < 256`.  The
audit source is pinned at SHA-256
`c48157f690d0f23fca66f267183e54d8d69a9d39853794dc824439936237baf8`.

## Source-equivalence boundary

The final and older replay manifests have the same 85 paths.  Seventy-seven
hashes agree.  The eight changed paths are exactly `Makefile`, `README.md`,
`REPLAY.md`, `clean_room_replay.py`, `paper.tex`,
`character_ring_iter/Makefile`,
`character_ring_iter/post_m29_bc_layered_mgf_mpfr.cpp`, and
`character_ring_iter/verify_chain.py`.  None of the H8--H24 frontier sources
changed.  Their seventeen C++ sources are separately authenticated by the
old and current `bc_interval_tail_code.sha256` manifests, agree byte-for-byte,
and retain byte-identical Makefile target recipes.  The current changed
layered-MGF source is covered by the final-source prefix replay; the current
chain verifier and hardened clean-driver preflight are covered by the
current-source stages and audits.

## Rebuild hardening

During the independent local rebuild, `make -B` exposed that `gen_header.py`
did not reproduce the manifested `lie_data.h`: it added unused B2, B3, and D4
namespaces.  No result from that build is accepted.  The generator now emits
exactly the checked-in namespaces and reproduces the pinned header SHA-256
`ba30755944b7afa31ff8d7de00ff1b3384274e3bacd59edadb5d3c6f1025f576`.
Both generator inputs are included in the character-pairing manifest.  The
distributed verifier regenerates the header in a temporary directory and
also checks that all 52 direct or generated file dependencies of the 47
`replay-build` targets are hash-manifested.

The initial eight-worker fleet launch is rejected for exceeding the intended
physical-memory envelope before any task completed.  Its preserved partial
logs are not inputs.  The accepted fleet candidate uses four workers and the
same 71 arithmetic task boundaries; details are in
`FULL_Q3_CLEAN_FLEET_REPLAY.md`.

The initial four-worker queue `fullq3cleanfleet0002` remains the authenticated
source of the 24 machine-A H25--H26 tasks.  To avoid making the proof audit
wait for its later tasks serially, the unchanged ledger is partitioned by
`run_full_q3_frontier_partition.py` (SHA-256
`784c47f797a61e805d509c8e6262ccf5abd7f3bf7757a3b8fb4dad6d84b80adf`).
Machine B receives seven H27 tasks; machine C receives the 40-task heavy
remainder.  The wrapper imports the pinned original ledger (SHA-256
`28b51756e053cb160e9cf5dde4f97a27ae416f084f7ef23ffedd02409098464a`)
and changes no C++ source, task boundary, or exact arithmetic.

Fleet build logs retain their original absolute temporary paths.  Archived
verification pins the digest and complete `source/character_ring_iter/...`
suffix of each build record and separately hashes the relocated archived file.
Thus relocation cannot weaken a source/binary identity check or make an
immutable valid build log spuriously unverifiable.

Machine A's build transcript also retains the historical Python-ledger digest
`5883f7d8961d34adde5ca78335480b8f412b30de35a4eadcd1a47bcfe3bb1922`.
That Python file is not an input to any recorded compiler command.  The
verifier pins this line as provenance, verifies every C++ source and resulting
binary independently, and requires the actual run top log, staged run ledger,
71-task specification, and every task command to use the authenticated run
ledger `28b51756...`.  The historical line is neither rewritten nor assigned a
proof role.

## Promotion condition

The fail-closed collector
`collect_full_q3_distributed_frontier.sh` (SHA-256
`7b147557743dc7c74c471d22a2c3da0741b552f800df312d6e1f09e32d901d5d`)
admits only the 24 names in machine A's frozen subset manifest and the complete
7- and 40-task B/C result directories.  It checks every per-partition digest,
script identity, current-source equality, terminal marker, and C resource
receipt; runs the full composition against hidden staging paths; then
atomically relocates the three partitions and source snapshot.  It next runs
the archived verifier against the archived source root and final fleet paths,
and finally reruns the live-root audit against those final paths.  Only after
all three passes does it write the top-level archive manifest.  The collector
itself pins the strengthened
verifier at SHA-256
`4097d664d30b50bb379aa1fdf1aef1ae930556c19214330283a75f5ae68f4aa4`
and the final audit wrapper at
`208062defb8d3686cb254950b88b48b3a5380e0e5d5ddecb33bd01d39ca22a24`.
It also archives and authenticates both monitor programs and their complete
watch logs.  Machine A's copied top transcript must hash to the digest written
in the freeze receipt, and every machine-C resource sample is parsed and
checked against the documented 12-GiB available-memory floor and 2-GiB
swap-use ceiling.  An unreviewed verifier replacement cannot promote the run.

The same atomic transaction archives an
`execution_source_snapshot/ginibre_q3` tree before the prepromotion
composition.  It reproduces and checks the 85-entry replay manifest, the
complete Part III source manifest, the accepted post-$m=29$ certificate
manifest, the current/legacy B/C code manifest, the character-pairing manifest,
the primary-reference manifest, and their recursive manifest closure.  The
collector also copies its own bytes, the composition verifier and wrapper,
both fleet drivers, the independent horizontal-strip oracle, and the
generated-header inputs.  The accepted archive reconstructs 601 files (about
34.5 MB), validates 19 embedded manifests and the snapshot-wide SHA-256
ledger, and runs the archived verifier through the complete source
reconstruction and full fleet composition.  The final collector stores that
successful full composition from the archived source root.  Thus a later
publication-only rebind of `paper.tex` cannot replace the source bytes against
which the 200-stage evidence was accepted.

The publication-only rebind is separately fail-closed by
`verify_full_q3_paper_rebind.py` (SHA-256
`fee2856f349b4896ae6495f8d51982c35acb1103ece111e87410200f2b63a145`,
also pinned by `certificates/full_q3/full_q3_paper_rebind_verifier.sha256`).
Its preflight pins the executed `paper.tex` digest
`6b2f273aa082f105df4bf60f8ba7047c62cafac88a3bc2c0a5cffdca95607d27`,
pairs all 850 numbered results with 850 proofs, and reconstructs the 2,099-edge
result-reference graph.  The baseline has exactly one nontrivial component:
`prop:post29`, `prop:post29-bc-half-bridge`, `prop:post29-direct`, and
`thm:classical`.  The verifier authorizes one exact prose replacement and no
other byte change.  Its in-memory candidate has digest
`d2fb8944ec081f47a1bfcf558bb65adb01936b0ac1d0dbbb1f561d63d7197010`,
removes only the edge `prop:post29 -> thm:classical`, retains 2,098 edges, and
has no nontrivial component.  The 348 results reachable from `thm:main` are
unchanged.  Mutated baseline bytes, an unchanged final paper, an extra final
byte, and a different edge edit are all rejected.
The same audit verifies every file in the 85-entry replay manifest and the
36-entry Part III manifest.  It permits only the induced `paper.tex` record
change in the former, and only the induced `paper.tex` and nested replay-
manifest record changes in the latter.  Their final required digests are
`bfe68eaf0ee4a5fe8263902b1145804a630381a54bfd4a9b4bb3135517196874`
and
`1cbd33b571e990c8b8bce47a1bd549129dea6cc89d965ec5951ec6f8f265fb59`.

The edit is intentionally not applied until the execution-source snapshot and
successful 200-stage transcript exist.  Afterwards the archived source root
is passed as `--execution-root` and the publication source root as
`--final-root`; the final PDFs, manifests, and dependency audit are then
rebuilt from the rebound source.

After the three fleet partitions have completed, the collection and exact
composition command is:

```text
ginibre_q3/collect_full_q3_distributed_frontier.sh
```

No distributed evidence is promoted until the composition transcript ends
with

```text
FULL_Q3_DISTRIBUTED theorem_scope=adjoint_generated all_positive_definite_obstruction=-8/279
FULL_Q3_DISTRIBUTED stages=200 prefix=49 legacy_h8_h24=34 fleet_a_h25_h26=24 fleet_b_h27_low=7 fleet_c_h27_h29_heavy=40 classical_post=30 exceptional=16 failures=0
FULL_Q3_DISTRIBUTED VERIFICATION: ALL PASS
```

The preceding audit-detail line must also contain
`frontier_certificate_files=7`, `frontier_certificate_blocks=71`, and
`certificate_data_projections=124`.

The exact paper-rebind audit, final document audit, complete manifest audit,
three-pass builds of both companion PDFs, Part III audit/build, and strict
referee disposition are run only after that marker is present and
authenticated.

## Publication finalization order

The following order is mandatory.  It prevents the publication-only edit from
changing the source snapshot against which the 200-stage computation was
accepted.

1. From `certificates/full_q3/distributed0001`, verify
   `distributed0001_manifest.sha256` with `sha256sum -c`.  Require the two
   composition transcripts and the collection transcript to retain their exact
   `ALL PASS` markers.
2. From `certificates/full_q3`, verify
   `full_q3_paper_rebind_verifier.sha256` with `sha256sum -c`.
3. Return to the `ginibre_q3` source root.  Before changing the publication
   tree, run

   ```text
   python3 verify_full_q3_paper_rebind.py \
     --execution-root certificates/full_q3/distributed0001/execution_source_snapshot/ginibre_q3 \
     --preflight
   ```

   and require `FULL_Q3_PAPER_REBIND VERIFICATION: ALL PASS`.
4. Apply only the fixed `OLD_PASSAGE` to `NEW_PASSAGE` replacement encoded in
   `verify_full_q3_paper_rebind.py`.  Change only the induced `paper.tex` entry
   in `replay_sources.sha256`, and only the induced `paper.tex` and nested
   replay-manifest entries in
   `certificates/full_q3/full_q3_source_manifest.sha256`.
5. Run the final-mode rebind audit:

   ```text
   python3 verify_full_q3_paper_rebind.py \
     --execution-root certificates/full_q3/distributed0001/execution_source_snapshot/ginibre_q3 \
     --final-root .
   ```

   It must report the fixed paper digest `d2fb8944...`, edge count
   `2099->2098`, component count `1->0`, unchanged reachability, and the two
   fixed final manifest digests.
6. Rerun the document, reference, caller-range, accepted-status, replay-input,
   and complete artifact audits from the final source.  Then build
   `paper.tex` and `paper_full.tex` for three clean `pdflatex` passes each and
   reject undefined references, undefined citations, duplicate labels, or a
   nonzero process status.
7. Run `make full-q3-extension` from `ginibre_q3`; require every arithmetic
   target, the 51-result document audit, the complete artifact audit, and the
   two-pass `full_q3_extension.pdf` build to succeed.
8. Verify that all three PDFs are nonempty and issue the unconditional referee
   disposition only after every preceding gate has passed.  Any failure reopens
   the exact failing obligation; it is not waived by the earlier conditional
   recommendation.

## Publication metadata, scope, and referee-correction rebind

After the preceding computation and dependency rebind were frozen, a
publication-only edit added named authorship, the source-audited discussion
of Sylvester (1980), Abdesselam (2022), and Chevyrev--Garban (2025), the
quantum-application boundary, and a code/data availability section.  It
changes no theorem, proof, arithmetic program, or certificate.  Its induced
`paper.tex`, replay-manifest, and Part III manifest changes are exact metadata
rebinds.  The final referee corrections additionally state the unresolved
central-positive-definite cone boundary, pin the audited computation by full
Git commit, make Parts I--II a formally inseparable submission component, add
the proposition-to-certificate map, and synchronize the publication date.
`verify_full_q3_related_work_rebind.py` reverses the fixed author
metadata, related-work and application-scope remarks, three bibliography
items, three source-table rows, the complete corrected availability section,
the central-cone sentence, both date changes, companion bibliography identity,
and all added document checks to recover the accepted pre-edit hashes.
The audit also verifies all
46 records in the final reference manifest and all 36 records in the final
Part III manifest.

From `certificates/full_q3`, first authenticate the new verifier:

```text
sha256sum -c full_q3_related_work_rebind_verifier.sha256
```

Then, from the `ginibre_q3` source root, run:

```text
python3 verify_full_q3_related_work_rebind.py
python3 verify_full_q3_document.py
sha256sum -c references/references_manifest.sha256
(cd certificates/full_q3 && sha256sum -c full_q3_source_manifest.sha256)
```

The related-work verifier must report the final Part III source digest
`ba089e58...0bf4e06`, 13 cited works, 46 reference records, 36 full-manifest
records, and `FULL_Q3_RELATED_WORK_REBIND VERIFICATION: ALL PASS`.  The final
Part III manifest digest is `6af35090...436a68`.  The earlier digest
`1cbd33b5...265fb59` remains the authenticated stage immediately after the
dependency rebind and is the baseline recovered by reversing this second,
publication-only edit.
