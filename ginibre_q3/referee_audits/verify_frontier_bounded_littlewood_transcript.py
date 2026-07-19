#!/usr/bin/env python3
"""Fail-closed audit of the full-scale bounded-Littlewood oracle transcript."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PRIME_RE = re.compile(
    r"^FRONTIER_DETERMINANT evaluated_prime=(\d+) prime_index=(\d+)/(\d+)$"
)
ROW_RE = re.compile(
    r"^B_DETERMINANT_BAD_COUNT h=(\d+) q=(\d+) rank=(\d+) j=(\d+) "
    r"finite=(\d+) stable=(\d+) bad=(\d+) accepted=(\d+) match=(\d+)$"
)
SUMMARY = (
    "FRONTIER_BOUNDED_LITTLEWOOD_ORACLE rows=33 cases=120 "
    "maximum_moment=121 primes=22 modulus_bits=682 failures=0"
)
PASS = "FRONTIER_BOUNDED_LITTLEWOOD_ORACLE VERIFICATION: ALL PASS"
EXIT = "__EXIT_STATUS=0"
EXPECTED_Q = {
    25: range(15, 44),
    26: range(15, 46),
    27: range(15, 48),
    28: range(15, 42),
}


def verify(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    prime_indices: dict[int, int] = {}
    rows: dict[tuple[int, int], tuple[int, ...]] = {}
    summaries = passes = exits = 0

    for number, line in enumerate(lines, 1):
        if match := PRIME_RE.fullmatch(line):
            prime, index, total = map(int, match.groups())
            if total != 22 or not 0 <= index < total:
                raise AssertionError(f"line {number}: invalid prime index")
            if index in prime_indices:
                raise AssertionError(f"line {number}: duplicate prime index {index}")
            if prime in prime_indices.values():
                raise AssertionError(f"line {number}: duplicate prime {prime}")
            prime_indices[index] = prime
            continue

        if match := ROW_RE.fullmatch(line):
            h, q, rank, j, finite, stable, bad, accepted, matched = map(
                int, match.groups()
            )
            key = (h, q)
            if key in rows:
                raise AssertionError(f"line {number}: duplicate row {key}")
            if rank != q - 1 or j != 2 * q + h:
                raise AssertionError(f"line {number}: wrong rank/moment coordinates")
            if not 0 <= finite <= stable or finite + bad != stable:
                raise AssertionError(f"line {number}: invalid exact subtraction")
            if bad != accepted or matched != 1:
                raise AssertionError(f"line {number}: accepted total mismatch")
            rows[key] = (rank, j, finite, stable, bad, accepted, matched)
            continue

        if line == SUMMARY:
            summaries += 1
        elif line == PASS:
            passes += 1
        elif line == EXIT:
            exits += 1
        else:
            raise AssertionError(f"line {number}: unknown transcript record: {line!r}")

    expected = {(h, q) for h, values in EXPECTED_Q.items() for q in values}
    if set(rows) != expected:
        missing = sorted(expected - set(rows))
        extra = sorted(set(rows) - expected)
        raise AssertionError(f"scope mismatch: missing={missing} extra={extra}")
    if set(prime_indices) != set(range(22)):
        raise AssertionError("prime-index coverage is not exactly 0..21")
    if (summaries, passes, exits) != (1, 1, 1):
        raise AssertionError(
            "terminal-record multiplicities are not exactly summary/pass/exit = 1/1/1"
        )
    maximum_stable_bits = max(values[3].bit_length() for values in rows.values())
    if maximum_stable_bits != 663:
        raise AssertionError(
            f"maximum stable bound changed: {maximum_stable_bits} bits, expected 663"
        )
    if 682 <= maximum_stable_bits:
        raise AssertionError("CRT modulus does not dominate the stable reconstruction bound")

    print(
        "FRONTIER_BOUNDED_LITTLEWOOD_TRANSCRIPT_AUDIT "
        f"rows={len(rows)} primes={len(prime_indices)} "
        f"maximum_stable_bits={maximum_stable_bits} failures=0"
    )
    print("FRONTIER_BOUNDED_LITTLEWOOD_TRANSCRIPT_AUDIT: ALL PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "transcript",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name(
            "frontier_bounded_littlewood_oracle_machine_c.log"
        ),
    )
    args = parser.parse_args()
    verify(args.transcript)


if __name__ == "__main__":
    main()
