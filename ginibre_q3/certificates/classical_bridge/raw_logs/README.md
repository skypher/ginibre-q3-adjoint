# Row-Gated Bridge Raw Logs

This directory vendors the raw logs cited by
`../row_gated_bridge_replay_transcript.log` for the m=26,27,28,29
row-gated bridge.

The subdirectory names preserve the historical `/tmp` basenames used in the
transcript:

- `ginibre_m24_logs`: copies of the two m24 reused-window inputs, recovered
  from the later `ginibre_m25_eval` cache.
- `ginibre_m25_eval`: reused-window inputs also cited by later rows.
- `ginibre_m26_logs`, `ginibre_m27_logs`, `ginibre_m28_logs`,
  `ginibre_m29_logs`: raw bridge logs for the four final row-gated slices.

`../raw_log_manifest.sha256` lists every vendored raw log.  It contains 297
entries and has SHA-256
`09b40fe1ebbdf8897471b513d6353c70eb1514238a76414cb20746e57d8310bd`.

The archive is partitioned by `../raw_log_classification.tsv`, with columns
`sha256`, `category`, `reason`, and `raw_logs/...` path.  Its SHA-256 is
`2277355c990e1e6c6b110bd98dcf9db47349e97d0bfe8c360e5441cdf2ab7e93`.
The category manifests are:

- `../accepted_log_manifest.sha256`: 124 status-zero replay, postcheck, or
  adjacent status-marker files used as row-acceptance evidence; SHA-256
  `975c4e143e1a9e273fea8eb074935860fcf26af9874be5da8bafdd33f4c4e6a5`.
- `../supporting_log_manifest.sha256`: 132 status-free arithmetic inputs or
  intermediate logs with no failure marker; SHA-256
  `21a8c98545124035c3abe42fd4ae1652f1255b217b4d3765dfac198402140155`.
- `../diagnostic_log_manifest.sha256`: 41 failed, interrupted, or superseded
  diagnostic logs containing a non-zero exit marker or explicit failure/missing
  marker; SHA-256
  `5007477325bf7ed035976d708b52417cb0e8cbf55191465dcc670b0f395cc34f`.

Only entries in the accepted manifest may serve as row-acceptance evidence.
Supporting entries may be cited as arithmetic inputs.  Diagnostic entries are
retained for audit trail only; they are not load-bearing certificate logs.
