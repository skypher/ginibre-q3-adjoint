#!/usr/bin/env python3
"""Assemble split D_b exact-boundary logs into Chain differences."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


LINE_RE = re.compile(
    r"D_(?P<rank>\d+) m_(?P<moment>\d+) \+= .*; Delta_\d+ = (?P<delta>-?\d+)"
)


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


def q3(moments: list[int], n: int) -> int:
    total = 0
    for k in range(n + 1):
        total += math.comb(n, k) * (
            moments[k + 2] * moments[n - k]
            - moments[k + 1] * moments[n - k + 1]
        )
    return 2 * total


def chain_diff(moments: list[int], chain_m: int) -> int:
    return q3(moments, 2 * chain_m + 3) - 4 * q3(moments, 2 * chain_m + 1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stable", type=Path, required=True)
    parser.add_argument("--chain-m", type=int, required=True)
    parser.add_argument("--rank-low", type=int, required=True)
    parser.add_argument("--rank-high", type=int, required=True)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()

    moment_max = 2 * args.chain_m + 5
    stable = read_stable(args.stable, moment_max)
    deltas: dict[tuple[int, int], int] = {}

    for log_path in args.logs:
        with log_path.open() as handle:
            for raw in handle:
                match = LINE_RE.search(raw)
                if not match:
                    continue
                rank = int(match.group("rank"))
                moment = int(match.group("moment"))
                if not (args.rank_low <= rank <= args.rank_high):
                    continue
                if moment > moment_max:
                    continue
                delta = int(match.group("delta"))
                key = (rank, moment)
                old = deltas.get(key)
                if old is not None and old != delta:
                    raise RuntimeError(
                        f"conflicting delta for D_{rank} m_{moment}: {old} vs {delta}"
                    )
                deltas[key] = delta

    for rank in range(args.rank_low, args.rank_high + 1):
        missing = [
            moment
            for moment in range(rank, moment_max + 1)
            if (rank, moment) not in deltas
        ]
        if missing:
            raise RuntimeError(f"D_{rank} missing moments: {missing[:12]}")
        row = stable[:]
        for moment in range(rank, moment_max + 1):
            row[moment] += deltas[(rank, moment)]
        diff = chain_diff(row, args.chain_m)
        print(f"D_{rank} exact assembled Chain diff D({args.chain_m}) = {diff}", flush=True)
        if diff <= 0:
            raise RuntimeError(f"D_{rank} assembled Chain diff is not positive")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
