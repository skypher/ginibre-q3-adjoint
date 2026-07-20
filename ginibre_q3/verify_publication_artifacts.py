#!/usr/bin/env python3
"""Authenticate the deterministic PDFs produced from the publication sources."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "publication_artifacts.sha256"
EXPECTED = {
    "paper.pdf",
    "full_q3_main.pdf",
    "full_q3_extension.pdf",
    "submission.pdf",
}
ROW = re.compile(r"([0-9a-f]{64})  ([^/]+\.pdf)")


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"PUBLICATION_ARTIFACTS: FAIL: {message}")


def main() -> int:
    if not MANIFEST.is_file():
        fail(f"missing {MANIFEST.name}")

    records: dict[str, str] = {}
    for number, line in enumerate(MANIFEST.read_text(encoding="ascii").splitlines(), 1):
        match = ROW.fullmatch(line)
        if match is None:
            fail(f"malformed manifest row {number}")
        expected_hash, name = match.groups()
        if name in records:
            fail(f"duplicate manifest entry {name}")
        records[name] = expected_hash

    if set(records) != EXPECTED:
        missing = sorted(EXPECTED - set(records))
        extra = sorted(set(records) - EXPECTED)
        fail(f"wrong manifest scope missing={missing} extra={extra}")

    for name in sorted(EXPECTED):
        source = ROOT / name
        if not source.is_file() or source.stat().st_size == 0:
            fail(f"missing or empty {name}")
        data = source.read_bytes()
        if not data.startswith(b"%PDF-"):
            fail(f"{name} is not a PDF")
        actual = hashlib.sha256(data).hexdigest()
        if actual != records[name]:
            fail(f"{name} sha256={actual} expected={records[name]}")
        print(f"PUBLICATION_ARTIFACTS file={name} bytes={len(data)} sha256={actual}")

    print("PUBLICATION_ARTIFACTS: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
