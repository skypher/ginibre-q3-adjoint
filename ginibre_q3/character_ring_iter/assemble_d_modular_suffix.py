#!/usr/bin/env python3
"""Assemble exact D_b suffix moments from modular row-capped replays."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


EXACT_RE = re.compile(
    r"D_(?P<rank>\d+) m_(?P<moment>\d+) \+= .*; Delta_\d+ = (?P<delta>-?\d+)"
)
MOD_RE = re.compile(
    r"D_(?P<rank>\d+) MOD (?P<prime>\d+) m_(?P<moment>\d+) "
    r"kept_even=(?P<kept>\d+) det=(?P<det>\d+)"
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


def chain_diff(moments: list[int], chain_m: int) -> int:
    total = 0
    for n_value, scale in ((2 * chain_m + 3, 1), (2 * chain_m + 1, -4)):
        for k in range(n_value + 1):
            c = 2 * scale * math.comb(n_value, k)
            total += c * (
                moments[k + 2] * moments[n_value - k]
                - moments[k + 1] * moments[n_value - k + 1]
            )
    return total


def crt_nonnegative(residues: list[tuple[int, int]], bound: int, label: str) -> int:
    value = 0
    modulus = 1
    for prime, residue in residues:
        if math.gcd(modulus, prime) != 1:
            raise RuntimeError(f"{label}: modulus factors are not coprime at {prime}")
        residue %= prime
        step = ((residue - value) % prime) * pow(modulus % prime, -1, prime)
        step %= prime
        value += modulus * step
        modulus *= prime
    if modulus <= bound:
        raise RuntimeError(
            f"{label}: CRT modulus product {modulus} does not exceed bound {bound}"
        )
    if value > bound:
        raise RuntimeError(
            f"{label}: reconstructed value {value} exceeds bound {bound}"
        )
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stable", type=Path, required=True)
    parser.add_argument("--chain-m", type=int, required=True)
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--exact-through", type=int, required=True)
    parser.add_argument("--suffix-low", type=int, required=True)
    parser.add_argument("--suffix-high", type=int, required=True)
    parser.add_argument("--prefix-log", type=Path, action="append", default=[])
    parser.add_argument("--mod-log", type=Path, action="append", required=True)
    parser.add_argument(
        "--det-only",
        action="store_true",
        help="reconstruct and print only determinant suffix deltas from modular logs",
    )
    args = parser.parse_args()

    moment_max = 2 * args.chain_m + 5
    stable = [] if args.det_only else read_stable(args.stable, moment_max)
    exact_delta: dict[int, int] = {}
    residues: dict[int, dict[int, tuple[int, int]]] = {}

    if not args.det_only:
        if not args.prefix_log:
            raise RuntimeError("non-det-only assembly requires at least one prefix log")
        for prefix_log in args.prefix_log:
            with prefix_log.open() as handle:
                for raw in handle:
                    match = EXACT_RE.search(raw)
                    if not match:
                        continue
                    rank = int(match.group("rank"))
                    moment = int(match.group("moment"))
                    if rank == args.rank and moment <= moment_max:
                        exact_delta[moment] = int(match.group("delta"))

    for mod_log in args.mod_log:
        with mod_log.open() as handle:
            for raw in handle:
                match = MOD_RE.search(raw)
                if not match:
                    continue
                rank = int(match.group("rank"))
                if rank != args.rank:
                    continue
                prime = int(match.group("prime"))
                moment = int(match.group("moment"))
                if args.suffix_low <= moment <= args.suffix_high:
                    residues.setdefault(moment, {})[prime] = (
                        int(match.group("kept")),
                        int(match.group("det")),
                    )

    missing_exact = [
        moment for moment in range(args.rank, args.exact_through + 1)
        if moment not in exact_delta
    ]
    if missing_exact:
        raise RuntimeError(f"D_{args.rank} missing exact prefix moments: {missing_exact[:12]}")

    missing_mod = [
        moment for moment in range(args.suffix_low, args.suffix_high + 1)
        if moment not in residues
    ]
    if missing_mod:
        raise RuntimeError(f"D_{args.rank} missing modular suffix moments: {missing_mod[:12]}")

    row_cap = 2 * args.rank
    path_branch_bound = math.comb(row_cap, 2)

    if args.det_only:
        for moment in range(args.suffix_low, args.suffix_high + 1):
            by_prime = residues[moment]
            det_residues = sorted((prime, pair[1]) for prime, pair in by_prime.items())
            det = crt_nonnegative(
                det_residues,
                path_branch_bound**moment,
                f"D_{args.rank} m_{moment} determinant",
            )
            print(f"D_{args.rank} Delta_{moment} = {det}", flush=True)
        return 0

    moments = stable[:]
    for moment in range(args.rank, args.exact_through + 1):
        moments[moment] = stable[moment] + exact_delta[moment]

    for moment in range(args.suffix_low, args.suffix_high + 1):
        by_prime = residues[moment]
        kept_residues = sorted((prime, pair[0]) for prime, pair in by_prime.items())
        det_residues = sorted((prime, pair[1]) for prime, pair in by_prime.items())
        if len(kept_residues) != len(det_residues):
            raise RuntimeError(f"D_{args.rank} m_{moment}: inconsistent residues")
        kept = crt_nonnegative(
            kept_residues,
            stable[moment],
            f"D_{args.rank} m_{moment} kept-even",
        )
        det = crt_nonnegative(
            det_residues,
            path_branch_bound**moment,
            f"D_{args.rank} m_{moment} determinant",
        )
        o_even = kept - stable[moment]
        delta = o_even + det
        moments[moment] = stable[moment] + delta
        print(f"D_{args.rank} O_even_{moment} = {o_even}", flush=True)
        print(f"D_{args.rank} Delta_{moment} = {det}", flush=True)
        print(
            f"  D_{args.rank} m_{moment} += O_even {o_even} + det {det}; "
            f"Delta_{moment} = {delta}; moment_{moment} = {moments[moment]}",
            flush=True,
        )

    diff = chain_diff(moments, args.chain_m)
    print(
        f"D_{args.rank} modular suffix exact Chain diff D({args.chain_m}) = {diff}",
        flush=True,
    )
    if diff <= 0:
        raise RuntimeError(f"D_{args.rank} modular suffix Chain diff is not positive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
