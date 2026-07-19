# Character-ring pairing source regeneration

Run id `ginibre_q3_character_pairing_20260713` independently regenerates the
load-bearing adjoint moments by exact Racah--Speiser iteration and character
orthogonality.  If

`chi_ad^j = sum_lambda c_j(lambda) chi_lambda`,

then exact character orthogonality gives

`m_(2j) = sum_lambda c_j(lambda)^2` and
`m_(2j+1) = sum_lambda c_j(lambda)c_(j+1)(lambda)`.

Consequently the E8 source through `m_100` requires iteration only through
`ad^50`.  The verifier also regenerates the full G2, F4, E6, and E7 moment
sources and every consumed classical B/C/D moment or stable-moment correction
by the same pairing identity.  It checks each Cartan matrix, reconstructs the
complete root set by simple reflections, verifies the adjoint root and
zero-weight multiplicities, rejects negative representation multiplicities,
and compares every exact GMP moment with the archived claim.  For a classical
correction, it also regenerates the stable sequence from the paper's proved
recurrence
`s_(n+1) = n s_n + n s_(n-1) - binom(n,2) s_(n-2)`, checks every entry in the
archived OEIS comparison table, and subtracts the regenerated rather than the
archived value.

## Primary E8 run

- Run id: `ginibre_q3_character_pairing_20260713_e8_c`
- Host: machine C, `nb1cb2f` (`149.50.101.50`)
- Assignment: E8, moments `m_0..m_100`, tensor powers `0..50`
- Threads: one; the exact map update is deterministic and memory bounded
- Remote work directory: `/tmp/ginibre_q3_character_pairing_20260713_e8_c`
- Accepted local log:
  `certificates/character_ring_pairing/e8_m0_m100_pairing_gmp_machine_c.log`
- Command:

  `./verify_character_ring_moment_sources_gmp E8 100 e8_70.log arxiv_2412_21189_e8_m71_m100.txt`

- Kill condition: any root-datum mismatch, negative character multiplicity,
  coordinate overflow, source mismatch, nonzero exit status, or memory
  exhaustion.  On a kill, no source table is promoted and the run switches to
  obstruction extraction.
- Promotion check: the accepted log must contain a summary beginning
  `SUMMARY group=E8 moments=101 source_values=101 failures=0`.

The other exceptional rows and all consumed classical source groups are short
enough for the ordinary clean-room replay on `optimus`; they are deliberately
regenerated there on every replay rather than authenticated from this
directory.  This machine-C log records an independent-hardware replication of
the heaviest single source-regeneration stage.

The adjacent manifest also pins `lie_data.py`, `gen_header.py`, and the exact
generated `lie_data.h`.  A forced regeneration must reproduce the manifested
header byte-for-byte before the pairing verifier is built.  The classical
`B`, `C`, and `D` rows are constructed algorithmically inside the verifier;
the generated header therefore contains only the exceptional data and the
small `A_1,A_2` audit fixtures that its checked-in version actually exposes.
