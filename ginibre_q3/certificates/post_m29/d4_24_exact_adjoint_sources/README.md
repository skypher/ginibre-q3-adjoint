# Exact type-D adjoint-moment sources

The top-level status-zero logs are the full signed finite-rank moments used by
the accepted D4–D24 bridge and terminal-tail replays.  Every active log was
generated on machine C (`nb1cb2f`) and records the SHA-256 of
`classical_boundary_certificate.cpp` before the generator output.  A source
row has the checked identity

`moment_j = stable_j + O_even_j + determinant_j`.

The `superseded_optimus/` subdirectory preserves earlier C++/GMP runs for
audit history only.  Active replay globs do not descend into it, and both
C++ consumers reject a full-moment source whose host is not machine C.

The bridge uses exact prefixes through m59 for D4 (independently reconstructed
from five modular Weyl runs), m44 for D5, m40 for D6–D7, m38 for D8–D18, m34
for D19, and m30 for D20–D24.  The D10/D11 terminal certificate uses full
adjoint moments through m38 and m46, respectively.  The D12–D24 exact-MGF
certificate uses the rank-dependent source degrees printed in its replay log.

The manifest in the parent directory authenticates every file.  The raw
modular Weyl inputs and their two independent D4 reconstruction summaries are
retained under `certificates/classical_bridge/raw_logs/ginibre_m27_logs/rank4`.

## Machine-C replacement run (2026-07-13)

Run id `d4_24_exact_moments_c_20260713` replaces the load-bearing source
shards that had been generated on `optimus`.  The primary host is machine C
(`nb1cb2f`, 32 hardware threads).  Existing completed shards cover
`D12..D17,m12..m38` and `D18..D24,m18..m30`; two independent processes cover
the remaining ranges `D18..D19,m31..m34` and `D18,m35..m36`.  Each process
uses the pinned C++/GMP generator source hash
`bc22efcafafe2d40b5f5430cb8f9cb046bc6b9b188da684e61dc6fdfaa27077e`
with `--reverse-target-sums --progress`.  Remote logs live under
`/tmp/ginibre_q3_dlow_root_20260713/logs/`.  The merge check requires exactly
one machine-C host line, one pinned source-hash line, status zero, all expected
rank/moment rows, and the exact identities
`O_even + determinant = Delta` and `stable + Delta = moment` before the parent
manifest is regenerated.

The same run id partitions the D11 terminal source into independent
machine-C processes for `m43`, `m44`, `m45`, and `m46` (one target moment per
process); completed machine-C logs for `m11..m42` are reused.  The merge check
is identical and the C++/GMP tail verifier additionally requires the complete,
gap-free interval `m11..m46` before proving its exact rational margin.
