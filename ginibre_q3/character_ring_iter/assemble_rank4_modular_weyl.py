#!/usr/bin/env python3
"""Assemble exact rank-four Weyl moments from modular constant-term replays."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import classical_tail_constants as tail_constants


MOD_RE = re.compile(
    r"(?P<family>[BCD])_4 MOD (?P<prime>\d+) m_(?P<moment>\d+) "
    r"moment=(?P<residue>\d+)"
)


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


def partial_chain_diff_lower(
    moments: list[int], stable: list[int], chain_m: int, partial_through: int
) -> tuple[int, list[tuple[int, int, int, int, str]]]:
    """Lower-bound D(chain_m) when later moments are only interval-bounded."""

    moment_max = 2 * chain_m + 5
    padded = moments + [0] * (moment_max - partial_through)
    base = chain_diff(padded, chain_m)
    lower = base
    choices: list[tuple[int, int, int, int, str]] = []
    for moment in range(partial_through + 1, moment_max + 1):
        probe = padded.copy()
        probe[moment] = 1
        coefficient = chain_diff(probe, chain_m) - base
        if coefficient < 0:
            lower += coefficient * stable[moment]
            choices.append(
                (moment, coefficient, stable[moment], coefficient * stable[moment], "S")
            )
        else:
            choices.append((moment, coefficient, stable[moment], 0, "0"))
    return lower, choices


def crt_nonnegative(residues: list[tuple[int, int]], bound: int, label: str) -> int:
    value = 0
    modulus = 1
    for prime, residue in residues:
        if math.gcd(modulus, prime) != 1:
            raise RuntimeError(f"{label}: modulus factors are not coprime at {prime}")
        residue %= prime
        step = ((residue - value) % prime) * pow(modulus % prime, -1, prime)
        value += modulus * (step % prime)
        modulus *= prime
    if modulus <= bound:
        raise RuntimeError(
            f"{label}: CRT modulus product {modulus} does not exceed bound {bound}"
        )
    if value > bound:
        raise RuntimeError(f"{label}: reconstructed value {value} exceeds bound {bound}")
    return value


def rank4_dimension(family: str) -> int:
    if family in {"B", "C"}:
        return 36
    if family == "D":
        return 28
    raise RuntimeError(f"unknown rank-four family {family}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("B", "C", "D"), required=True)
    parser.add_argument("--chain-m", type=int, required=True)
    parser.add_argument("--mod-log", type=Path, action="append", required=True)
    parser.add_argument(
        "--partial-through",
        type=int,
        help=(
            "reconstruct exact moments only through this index and use stable "
            "interval bounds for later moments"
        ),
    )
    args = parser.parse_args()

    moment_max = 2 * args.chain_m + 5
    reconstruct_through = args.partial_through
    if reconstruct_through is None:
        reconstruct_through = moment_max
    if reconstruct_through < 0 or reconstruct_through > moment_max:
        raise RuntimeError("--partial-through must lie in the moment range")
    if reconstruct_through < 2 * args.chain_m + 3:
        raise RuntimeError(
            "--partial-through must reach the top moment in Q3_M(2m+3)"
        )
    residues: dict[int, dict[int, int]] = {}

    for log_path in args.mod_log:
        with log_path.open() as handle:
            for raw in handle:
                match = MOD_RE.search(raw)
                if not match:
                    continue
                family = match.group("family")
                if family != args.family:
                    continue
                prime = int(match.group("prime"))
                moment = int(match.group("moment"))
                if 0 <= moment <= reconstruct_through:
                    residues.setdefault(moment, {})[prime] = int(match.group("residue"))

    missing = [
        moment for moment in range(reconstruct_through + 1)
        if moment not in residues
    ]
    if missing:
        raise RuntimeError(f"{args.family}_4 missing modular moments: {missing[:12]}")

    dimension = rank4_dimension(args.family)
    moments: list[int] = []
    for moment in range(reconstruct_through + 1):
        by_prime = sorted(residues[moment].items())
        value = crt_nonnegative(
            by_prime,
            dimension**moment,
            f"{args.family}_4 m_{moment}",
        )
        moments.append(value)
        if moment >= reconstruct_through - 12:
            print(f"{args.family}_4 m_{moment} = {value}", flush=True)

    if reconstruct_through < moment_max:
        stable = tail_constants.read_stable_moments(moment_max)
        lower, choices = partial_chain_diff_lower(
            moments, stable, args.chain_m, reconstruct_through
        )
        print(
            f"{args.family}_4 partial modular Weyl lower D({args.chain_m}) "
            f"through m_{reconstruct_through} = {lower}",
            flush=True,
        )
        for moment, coefficient, stable_value, contribution, endpoint in choices:
            print(
                f"{args.family}_4 interval moment m_{moment}: "
                f"coefficient={coefficient}; stable_bound={stable_value}; "
                f"endpoint={endpoint}; lower contribution={contribution}",
                flush=True,
            )
        if lower <= 0:
            raise RuntimeError(
                f"{args.family}_4 partial modular Weyl lower is not positive"
            )
        return 0

    diff = chain_diff(moments, args.chain_m)
    print(
        f"{args.family}_4 modular Weyl exact Chain diff D({args.chain_m}) = {diff}",
        flush=True,
    )
    if diff <= 0:
        raise RuntimeError(f"{args.family}_4 modular Weyl Chain diff is not positive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
