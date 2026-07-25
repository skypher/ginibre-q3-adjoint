#!/usr/bin/env python3
"""Authenticate the SU(2) research certificate packages.

These artifacts are not manifested publication sources, so nothing used to
check them.  Three separate truncation incidents reached the repository or
its branches before this existed:

  * a corrupted O11 frontier verifier source (deleted in ace5997),
  * a mismatched O13 package chunk (replaced in 4bc0bed),
  * a truncated O11 replay archive still sitting on
    agent/o11-replay-restore-20260725-v2, whose three chunks all stop at
    exactly the writer's per-call byte cap with no terminal partial chunk.

Every check here is on the *decoded* payload, so it is independent of
whether a package is stored as one file or split into chunks.  Splitting is
not required by GitHub in the first place: the largest payload is about 17
kB against a 50 MB warning threshold.

Fails closed two ways:
  1. a registered package that is missing, undecodable, or hash-mismatched;
  2. a package-shaped artifact in certificates/ that is not registered, so a
     new package cannot be added without recording its hash here.
"""

from __future__ import annotations

import base64
import binascii
import gzip
import hashlib
import re
import sys
from pathlib import Path

# name -> (stem, kind, sha256 of decoded payload)
#   kind "tar"    : base64 -> gzip tarball
#   kind "ledger" : base64 -> gzip -> Python literal
PACKAGES = {
    "o13_initial_frontier": (
        "su2_o13_initial_frontier_package.b64",
        "tar",
        "a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2",
    ),
    "o13_support245": (
        "su2_o13_support245_package.b64",
        "tar",
        "da31783b2297c692d7dbc4b10e9f41bc7a03b07b664b43c137403ad06995c76f",
    ),
    "odd_orbit_rank5_ledger": (
        "su2_odd_orbit_rank5_ledger",
        "ledger",
        "06709b4ef2a1fcecc9a2eddf52dd2028bad878fb5dafad4f7e48d44dc37df5ce",
    ),
}

# Artifacts that look like packages but are plain evidence, not payloads.
IGNORE = re.compile(r"\.log$|\.md$|\.tsv$|\.txt$|\.sha256$|\.csv$")


class PackageFailure(RuntimeError):
    pass


def read_payload(certdir: Path, stem: str) -> bytes:
    """Concatenate chunks if present, else read the single file.

    Never mixes the two layouts, and never hardcodes a chunk count -- the
    previous rank-five reader hardcoded range(4), which would have silently
    dropped a fifth chunk.
    """
    chunks = sorted(certdir.glob(f"{stem}.part*"))
    if chunks:
        expected = [f"{stem}.part{i:02d}" for i in range(len(chunks))]
        actual = [c.name for c in chunks]
        if actual != expected:
            raise PackageFailure(
                f"{stem}: non-contiguous chunks {actual}, expected {expected}"
            )
        return b"".join(c.read_bytes() for c in chunks)
    single = certdir / stem
    if single.is_file():
        return single.read_bytes()
    raise PackageFailure(f"{stem}: no chunks and no single file")


def decode(stem: str, kind: str, payload: bytes) -> bytes:
    # Chunk files carry line breaks; strip all whitespace before validating so
    # that a genuinely non-base64 byte still fails rather than being ignored.
    compact = b"".join(payload.split())
    try:
        raw = base64.b64decode(compact, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise PackageFailure(f"{stem}: base64 decode failed ({exc})") from exc
    if kind == "tar":
        try:
            gzip.decompress(raw)
        except (OSError, EOFError) as exc:
            raise PackageFailure(
                f"{stem}: payload is not a complete gzip stream ({exc}); "
                "this is the truncation signature"
            ) from exc
    elif kind == "ledger":
        try:
            gzip.decompress(raw).decode()
        except (OSError, EOFError, UnicodeDecodeError) as exc:
            raise PackageFailure(f"{stem}: ledger did not decompress ({exc})") from exc
    return raw


def main() -> int:
    certdir = Path(__file__).resolve().parent / "certificates"
    failures: list[str] = []
    checked = 0

    for name, (stem, kind, want) in sorted(PACKAGES.items()):
        try:
            payload = read_payload(certdir, stem)
            raw = decode(stem, kind, payload)
        except PackageFailure as exc:
            failures.append(str(exc))
            continue
        got = hashlib.sha256(raw).hexdigest()
        if got != want:
            failures.append(f"{stem}: sha256 {got}, expected {want}")
            continue
        checked += 1
        print(f"SU2_PACKAGE {name} bytes={len(raw)} sha256={got[:16]}... OK")

    # Fail closed on unregistered package-shaped artifacts.
    known: set[str] = set()
    for stem, _, _ in PACKAGES.values():
        known.add(stem)
        known.update(p.name for p in certdir.glob(f"{stem}.part*"))
    for path in sorted(certdir.iterdir()):
        if not path.is_file() or IGNORE.search(path.name):
            continue
        if ".b64" in path.name or ".part" in path.name:
            if path.name not in known:
                failures.append(
                    f"{path.name}: package-shaped artifact is not registered "
                    "in verify_su2_certificate_packages.py"
                )

    if failures:
        print()
        for f in failures:
            print(f"SU2_PACKAGE FAILURE: {f}", file=sys.stderr)
        print(
            f"\nSU2_CERTIFICATE_PACKAGES VERIFICATION: FAILED "
            f"({len(failures)} problem(s))",
            file=sys.stderr,
        )
        return 1

    print(f"SU2_CERTIFICATE_PACKAGES packages={checked}")
    print("SU2_CERTIFICATE_PACKAGES VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
