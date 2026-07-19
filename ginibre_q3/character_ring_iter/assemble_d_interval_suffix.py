#!/usr/bin/env python3
"""Assemble D_b prefix-exact plus determinant suffix interval Chain bounds."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


EXACT_RE = re.compile(
    r"D_(?P<rank>\d+) m_(?P<moment>\d+) \+= .*; Delta_\d+ = (?P<delta>-?\d+)"
)
DET_RE = re.compile(r"D_(?P<rank>\d+) Delta_(?P<moment>\d+) = (?P<delta>-?\d+)")
OEVEN_RE = re.compile(r"D_(?P<rank>\d+) O_even_(?P<moment>\d+) = (?P<delta>-?\d+)")


def read_stable(path: Path, moment_max: int) -> list[int]:
    values: list[int | None] = [None] * (moment_max + 1)
    with path.open() as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) < 2:
                continue
            index = int(fields[0])
            if 0 <= index <= moment_max:
                values[index] = int(fields[1])
    missing = [index for index, value in enumerate(values) if value is None]
    if missing:
        raise RuntimeError(f"missing stable moments: {missing}")
    return [int(value) for value in values]


def interval_product_lower(
    lower: list[int],
    upper: list[int],
    first: int,
    second: int,
    coefficient: int,
) -> int:
    products = (
        lower[first] * lower[second],
        lower[first] * upper[second],
        upper[first] * lower[second],
        upper[first] * upper[second],
    )
    if coefficient >= 0:
        return coefficient * min(products)
    return coefficient * max(products)


def chain_interval_lower(lower: list[int], upper: list[int], chain_m: int) -> int:
    total = 0
    for n_value, scale in ((2 * chain_m + 3, 1), (2 * chain_m + 1, -4)):
        for k in range(n_value + 1):
            c = 2 * scale * math.comb(n_value, k)
            total += interval_product_lower(lower, upper, k + 2, n_value - k, c)
            total += interval_product_lower(lower, upper, k + 1, n_value - k + 1, -c)
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stable", type=Path, required=True)
    parser.add_argument("--chain-m", type=int, required=True)
    parser.add_argument("--rank-low", type=int, required=True)
    parser.add_argument("--rank-high", type=int, required=True)
    parser.add_argument("--exact-through", type=int, required=True)
    parser.add_argument("--suffix-low", type=int, required=True)
    parser.add_argument("--suffix-high", type=int, required=True)
    parser.add_argument("--prefix-log", type=Path, action="append", required=True)
    parser.add_argument(
        "--allow-nonpositive",
        action="store_true",
        help="Diagnostic mode: print every interval lower bound without failing at the first nonpositive row.",
    )
    parser.add_argument(
        "--diagnose-collapses",
        action="store_true",
        help="For each row, also print the lower bound after collapsing each suffix moment to either endpoint.",
    )
    parser.add_argument("det_logs", type=Path, nargs="*")
    args = parser.parse_args()

    moment_max = 2 * args.chain_m + 5
    stable = read_stable(args.stable, moment_max)
    exact: dict[tuple[int, int], int] = {}
    det: dict[tuple[int, int], int] = {}
    oeven: dict[tuple[int, int], int] = {}

    for prefix_log in args.prefix_log:
        with prefix_log.open() as handle:
            for raw in handle:
                match = EXACT_RE.search(raw)
                if match:
                    rank = int(match.group("rank"))
                    moment = int(match.group("moment"))
                    if args.rank_low <= rank <= args.rank_high and moment <= moment_max:
                        exact[(rank, moment)] = int(match.group("delta"))
                        continue
                match = OEVEN_RE.search(raw)
                if match:
                    rank = int(match.group("rank"))
                    moment = int(match.group("moment"))
                    if args.rank_low <= rank <= args.rank_high and moment <= moment_max:
                        oeven[(rank, moment)] = int(match.group("delta"))
                    continue

    for log_path in args.det_logs:
        with log_path.open() as handle:
            for raw in handle:
                match = DET_RE.search(raw)
                if not match:
                    continue
                rank = int(match.group("rank"))
                moment = int(match.group("moment"))
                if args.rank_low <= rank <= args.rank_high:
                    det[(rank, moment)] = int(match.group("delta"))

    for rank in range(args.rank_low, args.rank_high + 1):
        missing_exact = [
            moment for moment in range(rank, args.exact_through + 1)
            if (rank, moment) not in exact
        ]
        missing_det = []
        missing_det = [
            moment for moment in range(args.suffix_low, args.suffix_high + 1)
            if (rank, moment) not in det and (rank, moment) not in exact
        ]
        if missing_exact:
            raise RuntimeError(f"D_{rank} missing exact prefix moments: {missing_exact[:12]}")
        if missing_det:
            raise RuntimeError(f"D_{rank} missing determinant suffix moments: {missing_det[:12]}")

        lower = stable[:]
        upper = stable[:]
        for moment in range(rank, args.exact_through + 1):
            value = stable[moment] + exact[(rank, moment)]
            lower[moment] = value
            upper[moment] = value
        for moment in range(args.suffix_low, args.suffix_high + 1):
            if (rank, moment) in exact:
                value = stable[moment] + exact[(rank, moment)]
                lower[moment] = value
                upper[moment] = value
            elif (rank, moment) in oeven:
                value = stable[moment] + oeven[(rank, moment)] + det[(rank, moment)]
                lower[moment] = value
                upper[moment] = value
            else:
                det_value = det[(rank, moment)]
                lower[moment] = det_value
                upper[moment] = stable[moment] + det_value

        bound = chain_interval_lower(lower, upper, args.chain_m)
        print(
            f"D_{rank} suffix interval lower Chain diff D({args.chain_m}) = {bound}",
            flush=True,
        )
        if args.diagnose_collapses:
            candidates: list[tuple[int, int, str, int]] = []
            for moment in range(args.suffix_low, args.suffix_high + 1):
                for endpoint_name, endpoint_value in (
                    ("det", lower[moment]),
                    ("stable+det", upper[moment]),
                ):
                    trial_lower = lower[:]
                    trial_upper = upper[:]
                    trial_lower[moment] = endpoint_value
                    trial_upper[moment] = endpoint_value
                    trial_bound = chain_interval_lower(trial_lower, trial_upper, args.chain_m)
                    candidates.append((trial_bound - bound, moment, endpoint_name, trial_bound))
            candidates.sort(reverse=True)
            for improvement, moment, endpoint_name, trial_bound in candidates:
                print(
                    f"  collapse m_{moment} to {endpoint_name}: "
                    f"bound {trial_bound} (improvement {improvement})",
                    flush=True,
                )
        if bound <= 0 and not args.allow_nonpositive:
            raise RuntimeError(f"D_{rank} interval lower is not positive")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
