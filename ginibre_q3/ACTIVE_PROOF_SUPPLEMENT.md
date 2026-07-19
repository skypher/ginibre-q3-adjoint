# Active proof supplement

This file is the minimal index of the results and computations used by the
publication theorem. It is an index, not a substitute for the numbered
proofs. A diagnostic, conditional, exploratory, or superseded result in
`paper_full.pdf` is outside the proof unless it appears below as an explicit
incoming node.

## Submission components

| Component | Active role |
|---|---|
| `paper.pdf` | Parts I--II theorem, reductions, compact family closures, and final assembly |
| `paper_full.pdf` | The 58 B/C correction-prefix proofs and the half-stable bridge consumed by the compact residual closure |
| `full_q3_extension.pdf` | Part III cone reduction, higher-minus hierarchy, finite certificates, and final cone promotion |

## Parts I--II active chain

1. `lem:moment-formula`, `lem:even`, and `lem:chain-propagation` reduce the
   two-minus integral to adjoint moments and Chain inequalities.
2. `thm:type-a` consumes the stable-rank theorem, the low-rank ratio
   theorems, and their exact overlap certificates.
3. `prop:bc-active-correction-prefix-contract` identifies exactly 58 active
   B/C correction results: offsets `0..28` for each family. The compact
   residual proof consumes offsets `0..27`; offset `28` is overlap.
4. `prop:post29-bc-half-bridge` converts the active correction intervals into
   the common stable-moment lower bound used by
   `prop:post29-bc-residual-closure`.
5. `thm:classical` consumes that B/C closure, the type-D exact bridges, and
   the exact or outward-rounded tail suppliers.
6. `thm:exceptional` consumes the five exact character-ring prefixes and
   their bridge/tail certificates.
7. `lem:global-form`, classification, and the family theorems assemble
   `thm:main`.

The machine-checked caller contract is `verify_parts_i_ii_contract.py`. It
requires 58 active correction references, maximum consumed offset 27,
overlap offset 28, and the named half-bridge input. The complete arithmetic
regeneration command is `make -C ginibre_q3 clean-room-replay`.

## Part III active chain

| Result | Incoming mathematical or computational nodes |
|---|---|
| `thm:full-cone-reduction` | Generator hierarchy, multiplicative/convex closure, and uniform closure |
| `thm:two-minus-import` | `thm:main` from Parts I--II |
| `thm:type-a-full-cone` | Stable type-A law, exact SU(2)--SU(5) prefixes, and the SU(N≥6) residual certificate |
| `thm:exceptional-full-cone` | Exact exceptional moments, polynomial tails, and 1,265 residual checks |
| `thm:stable-bcd-full` | Stable Gaussian-square law and exact cumulant recurrence |
| `thm:bcd-full-hierarchy` | High/middle outward-MPFR suppliers; 58 low-tail schedules; 12,993 bounded-Littlewood residual checks; 4,869 B/D residual checks |
| `thm:full-adjoint-generated-q3` | Classification plus `thm:full-cone-reduction` |

The Part III arithmetic regeneration command is
`make -C ginibre_q3 full-q3-extension`; it includes the compressed independent
target `full-q3-bcd-independent-audit`. That target keeps the structurally
independent reverse-Pieri comparison through degree 40 (2,845 moments), then
checks the 538 required higher moments modulo 21 independent primes by a
plain-arithmetic determinant implementation using Newton differences. The
628-bit modulus exceeds the 607-bit character bound, so the residues determine
those moments uniquely. It finally checks the remaining 5,509 hierarchy
values exactly.

## Trusted-computation boundary

All sign decisions use GMP exact arithmetic or outward-rounded MPFR
intervals. Exact arithmetic removes rounding ambiguity but does not verify an
implementation. The compressed checker is separate source from the promotion
verifier, uses a different interpolation and modular-arithmetic path, and
fails closed on the moment, modulus-bound, overlap, or hierarchy scope.

`publication-preflight` checks structure, formulas, coverage, and manifests;
it does not recompute all signs. A complete final-source reproduction requires
both expensive aggregates. `run_final_publication_replay.sh` runs those
aggregates from a clean commit, records environment and source identities,
and rejects any replay that changes the source tree.

## Explicit exclusions

The following are not theorem premises: the author self-audit,
`referee_audits/`, historical route diaries, launch-failure logs, superseded
transcripts, the false `DEnv` route, offset-29 conditional continuations, and
any certificate not classified as accepted by the replay contract.
