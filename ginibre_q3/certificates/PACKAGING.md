# Packaging replay archives

## Do not split packages into chunks

Every SU(2) replay package in this directory is a few kilobytes. The largest
decoded payload is about 17 kB. GitHub warns at 50 MB and refuses at 100 MB,
so these are three orders of magnitude below any limit that would justify
splitting.

The chunk boundaries that have appeared here — 12001, 5000, and 4001 bytes —
are not GitHub limits. They are per-call write limits of the tool uploading
the file. Splitting to satisfy a writer's limit is what has corrupted these
archives four times:

| Incident | Outcome |
|---|---|
| O11 frontier verifier source | binary garbage mid-literal, deleted in `ace5997` |
| O13 package chunk | mismatched, replaced in `4bc0bed` |
| O11 replay archive | truncated at the write cap; still broken on `agent/o11-replay-restore-20260725-v2` |
| O13 support 259 package | never landed; transcript has no reproduction path |

The failure is always the same shape: a write stops at the cap, no terminal
partial chunk is produced, and nothing notices because the concatenation still
looks plausible. Store the package as **one file**.

## Verify before and after upload

`ginibre_q3/verify_su2_certificate_packages.py` authenticates every registered
package by decoding it and comparing a SHA-256 of the *decoded* payload, so it
is independent of storage layout. It runs on every push via
`.github/workflows/preflight.yml`.

It fails closed in both directions: a registered package that is missing,
undecodable, or hash-mismatched fails, and a package-shaped file that is not
registered also fails. So adding a new package requires adding its hash to
`PACKAGES` in that script — which is the point. A package nobody registered is
a package nobody checked.

## Adding a package

1. Build the archive and record `sha256sum` of the **decoded** payload.
2. Commit it as a single `.b64` file.
3. Add an entry to `PACKAGES` in `verify_su2_certificate_packages.py`.
4. Run `python3 ginibre_q3/verify_su2_certificate_packages.py` locally.

If a package genuinely will not fit in one write, that is a tooling problem to
fix at the writer, not by splitting the artifact. Readers still accept the
chunked layout for backward compatibility, but nothing should produce it.
