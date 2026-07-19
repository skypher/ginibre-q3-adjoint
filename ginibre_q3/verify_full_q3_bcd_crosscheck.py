#!/usr/bin/env python3
"""Compare independent reverse-Pieri and bounded-Littlewood moment logs."""

from __future__ import annotations

import argparse
from pathlib import Path
import re


REVERSE = re.compile(r"^FULL_Q3_MOMENT row=([BCD]_\d+) degree=(\d+) value=(\d+)$")
BOUNDED = re.compile(r"^FULL_Q3_BOUNDED_MOMENT row=([BCD]_\d+) degree=(\d+) value=(\d+)$")


def read_moments(path: Path, pattern: re.Pattern[str], maximum: int) -> dict[tuple[str, int], int]:
    values: dict[tuple[str, int], int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.fullmatch(line.strip())
        if match is None:
            continue
        degree = int(match.group(2))
        if degree > maximum:
            continue
        key = (match.group(1), degree)
        if key in values:
            raise SystemExit(f"FULL_Q3_BCD_CROSSCHECK: FAIL: duplicate {key} in {path}")
        values[key] = int(match.group(3))
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reverse_log", type=Path)
    parser.add_argument("bounded_log", type=Path)
    parser.add_argument("--maximum", type=int, default=59)
    args = parser.parse_args()
    if args.maximum < 7 or args.maximum > 59:
        raise SystemExit("FULL_Q3_BCD_CROSSCHECK: FAIL: maximum must lie in 7..59")

    reverse = read_moments(args.reverse_log, REVERSE, args.maximum)
    bounded = read_moments(args.bounded_log, BOUNDED, args.maximum)
    if not reverse or not bounded:
        raise SystemExit("FULL_Q3_BCD_CROSSCHECK: FAIL: an input contains no moments")
    if reverse.keys() != bounded.keys():
        missing_reverse = sorted(bounded.keys() - reverse.keys())[:5]
        missing_bounded = sorted(reverse.keys() - bounded.keys())[:5]
        raise SystemExit(
            "FULL_Q3_BCD_CROSSCHECK: FAIL: scope mismatch "
            f"missing_reverse={missing_reverse} missing_bounded={missing_bounded}"
        )
    mismatches = [key for key in reverse if reverse[key] != bounded[key]]
    if mismatches:
        key = mismatches[0]
        raise SystemExit(
            "FULL_Q3_BCD_CROSSCHECK: FAIL: value mismatch "
            f"key={key} reverse={reverse[key]} bounded={bounded[key]}"
        )
    maximum_seen = max(degree for _row, degree in reverse)
    if maximum_seen != args.maximum:
        raise SystemExit(
            "FULL_Q3_BCD_CROSSCHECK: FAIL: requested endpoint absent "
            f"requested={args.maximum} observed={maximum_seen}"
        )
    print(
        "FULL_Q3_BCD_CROSSCHECK "
        f"moments={len(reverse)} maximum={maximum_seen} implementations=2"
    )
    print("FULL_Q3_BCD_CROSSCHECK VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
