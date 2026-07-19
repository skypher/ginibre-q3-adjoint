#!/usr/bin/env python3
"""Exact B/C/D adjoint Chain bridge replay.

This script uses the local Racah-Speiser/Klimyk character-ring iterator to
compute adjoint moments and verify

    Q3(2m+3) - 4 Q3(2m+1) >= 0

for the nonstable classical ranks left by the stable-rank theorem in
CLASSICAL_M3_ONE_PROOF.md, Corollary 16.
"""

from __future__ import annotations

import argparse
import sys
import time
from math import comb
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import lie_data  # noqa: E402
from character_ring import compute_moments  # noqa: E402


def cartan_Bn(n: int) -> list[list[int]]:
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for i in range(n - 1):
        if i < n - 2:
            A[i][i + 1] = -1
            A[i + 1][i] = -1
        else:
            A[i][i + 1] = -2
            A[i + 1][i] = -1
    return A


def cartan_Cn(n: int) -> list[list[int]]:
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for i in range(n - 1):
        if i < n - 2:
            A[i][i + 1] = -1
            A[i + 1][i] = -1
        else:
            A[i][i + 1] = -1
            A[i + 1][i] = -2
    return A


def cartan_Dn(n: int) -> list[list[int]]:
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    A[0][2] = -1
    A[2][0] = -1
    A[1][2] = -1
    A[2][1] = -1
    for i in range(2, n - 1):
        A[i][i + 1] = -1
        A[i + 1][i] = -1
    return A


def install_classical_cartan(family: str, rank: int) -> str:
    group = f"{family}{rank}"
    if family == "B":
        lie_data.CARTAN_MATRICES[group] = cartan_Bn(rank)
    elif family == "C":
        lie_data.CARTAN_MATRICES[group] = cartan_Cn(rank)
    elif family == "D":
        if rank < 4:
            raise ValueError("D rank must be at least 4")
        lie_data.CARTAN_MATRICES[group] = cartan_Dn(rank)
    else:
        raise ValueError(f"unknown family {family}")
    return group


def q3(moments: list[int], n: int) -> int:
    return 2 * sum(
        comb(n, k)
        * (moments[k + 2] * moments[n - k] - moments[k + 1] * moments[n - k + 1])
        for k in range(n + 1)
    )


def chain_diff(moments: list[int], m: int) -> int:
    return q3(moments, 2 * m + 3) - 4 * q3(moments, 2 * m + 1)


def nonstable_ranks(family: str, max_chain: int) -> range:
    if family == "B":
        return range(2, max_chain + 2)
    if family == "C":
        return range(2, max_chain + 2)
    if family == "D":
        return range(4, 2 * max_chain + 6)
    raise ValueError(f"unknown family {family}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-chain", type=int, default=4)
    parser.add_argument("--families", default="B,C,D")
    parser.add_argument(
        "--rank-cap",
        type=int,
        default=None,
        help="optional cap for quick partial runs",
    )
    parser.add_argument(
        "--rank-min",
        type=int,
        default=None,
        help="optional lower rank bound for resumable row checks",
    )
    parser.add_argument(
        "--moment-progress",
        action="store_true",
        help="print each computed adjoint moment inside every rank replay",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    families = [x.strip() for x in args.families.split(",") if x.strip()]
    moment_max = 2 * args.max_chain + 5
    print(
        f"Exact classical Chain bridge: max_chain={args.max_chain}, "
        f"moment_max={moment_max}, families={','.join(families)}, "
        f"rank_min={args.rank_min}, rank_cap={args.rank_cap}",
        flush=True,
    )

    all_ok = True
    rows = []
    for family in families:
        ranks = list(nonstable_ranks(family, args.max_chain))
        if args.rank_min is not None:
            ranks = [rank for rank in ranks if rank >= args.rank_min]
        if args.rank_cap is not None:
            ranks = [rank for rank in ranks if rank <= args.rank_cap]
        for rank in ranks:
            group = install_classical_cartan(family, rank)
            t0 = time.time()
            moments = compute_moments(group, moment_max, verbose=args.moment_progress)
            diffs = [chain_diff(moments, m) for m in range(args.max_chain + 1)]
            q_values = [q3(moments, 2 * m + 1) for m in range(args.max_chain + 2)]
            ok = all(diff >= 0 for diff in diffs) and all(q >= 0 for q in q_values)
            all_ok = all_ok and ok
            rows.append((group, rank, min(diffs), diffs[-1], q_values[-1], time.time() - t0))
            print(
                f"{group}: min_diff_m0..{args.max_chain}={min(diffs)} "
                f"diff_m{args.max_chain}={diffs[-1]} "
                f"Q3({2 * args.max_chain + 3})={q_values[-1]} "
                f"m0..6={moments[:7]} "
                f"dt={time.time() - t0:.3f}s {'OK' if ok else 'FAIL'}",
                flush=True,
            )

    print("\nSummary table:", flush=True)
    print("group,rank,min_chain_diff,last_chain_diff,last_odd_Q3,seconds", flush=True)
    for group, rank, min_diff, last_diff, last_q3, seconds in rows:
        print(f"{group},{rank},{min_diff},{last_diff},{last_q3},{seconds:.3f}", flush=True)
    print(f"\nRESULT: {'ALL PASS' if all_ok else 'FAIL'}", flush=True)
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
