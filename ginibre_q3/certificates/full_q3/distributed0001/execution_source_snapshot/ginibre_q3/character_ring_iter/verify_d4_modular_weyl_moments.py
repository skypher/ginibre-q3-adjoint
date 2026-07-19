#!/usr/bin/env python3
"""Reconstruct the exact D4 adjoint moments from archived modular Weyl runs."""

from __future__ import annotations

import argparse
import hashlib
import math
import re
from pathlib import Path


MODULAR_ROW = re.compile(
    r"^D_4 MOD (\d+) m_(\d+) moment=(\d+)$"
)
EXACT_ROW = re.compile(r"^D_4 m_(\d+) = (\d+)$")
SOURCE_MOMENT_ROW = re.compile(
    r"D_4 m_(\d+) .*; moment_\d+ = (\d+)$"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_prime_u64(value: int) -> bool:
    if value < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small_primes:
        if value % prime == 0:
            return value == prime
    odd_part = value - 1
    power_of_two = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        power_of_two += 1
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(power_of_two - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def parse_modular_log(path: Path, maximum: int) -> tuple[int, list[int]]:
    rows: dict[int, int] = {}
    modulus: int | None = None
    for line in path.read_text(errors="strict").splitlines():
        match = MODULAR_ROW.fullmatch(line.strip())
        if match is None:
            continue
        printed_modulus, index, residue = map(int, match.groups())
        if modulus is None:
            modulus = printed_modulus
        if printed_modulus != modulus:
            raise SystemExit(f"mixed moduli in {path}")
        if index in rows and rows[index] != residue:
            raise SystemExit(f"conflicting residue m_{index} in {path}")
        if residue < 0 or residue >= printed_modulus:
            raise SystemExit(f"out-of-range residue m_{index} in {path}")
        rows[index] = residue
    if modulus is None or not is_prime_u64(modulus):
        raise SystemExit(f"missing or nonprime modulus in {path}")
    expected = set(range(maximum + 1))
    if set(rows) != expected:
        missing = sorted(expected - set(rows))
        extra = sorted(set(rows) - expected)
        raise SystemExit(f"incomplete modular moment log {path}: missing={missing} extra={extra}")
    return modulus, [rows[index] for index in range(maximum + 1)]


def parse_exact_summary(path: Path, first: int, maximum: int) -> dict[int, int]:
    rows: dict[int, int] = {}
    for line in path.read_text(errors="strict").splitlines():
        match = EXACT_ROW.fullmatch(line.strip())
        if match is not None:
            index, value = map(int, match.groups())
            rows[index] = value
    expected = set(range(first, maximum + 1))
    if set(rows) != expected:
        raise SystemExit(f"incomplete exact reconstruction summary: {path}")
    return rows


def parse_source_cross_check(paths: list[Path], maximum: int) -> dict[int, int]:
    rows: dict[int, int] = {}
    for path in paths:
        for line in path.read_text(errors="strict").splitlines():
            match = SOURCE_MOMENT_ROW.search(line.strip())
            if match is None:
                continue
            index, value = map(int, match.groups())
            if index > maximum:
                continue
            if index in rows and rows[index] != value:
                raise SystemExit(f"conflicting D4 source moment m_{index}")
            rows[index] = value
    return rows


def read_stable(path: Path, maximum: int) -> list[int]:
    rows: dict[int, int] = {}
    for line in path.read_text(errors="strict").splitlines():
        if not line.strip():
            continue
        fields = line.split()
        index = int(fields[0])
        if index <= maximum:
            rows[index] = int(fields[1])
    if set(rows) != set(range(maximum + 1)):
        raise SystemExit("stable moment table is incomplete")
    return [rows[index] for index in range(maximum + 1)]


def crt(residues: list[int], moduli: list[int]) -> int:
    value = 0
    modulus_product = 1
    for residue, modulus in zip(residues, moduli, strict=True):
        multiplier = ((residue - value) * pow(modulus_product, -1, modulus)) % modulus
        value += modulus_product * multiplier
        modulus_product *= modulus
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mod-log", action="append", type=Path, required=True)
    parser.add_argument("--reconstruction-log", type=Path, required=True)
    parser.add_argument("--reconstruction-check-log", type=Path, required=True)
    parser.add_argument("--source-cross-check-log", action="append", type=Path, default=[])
    parser.add_argument("--stable-moments", type=Path, required=True)
    parser.add_argument("--maximum", type=int, default=59)
    args = parser.parse_args()

    if len(args.mod_log) < 2 or args.maximum < 0:
        raise SystemExit("at least two modular logs and a nonnegative maximum are required")
    parsed = [parse_modular_log(path, args.maximum) for path in args.mod_log]
    moduli = [item[0] for item in parsed]
    if len(set(moduli)) != len(moduli):
        raise SystemExit("duplicate CRT modulus")
    if any(math.gcd(first, second) != 1 for i, first in enumerate(moduli)
           for second in moduli[i + 1:]):
        raise SystemExit("CRT moduli are not pairwise coprime")

    modulus_product = math.prod(moduli)
    dimension_bound = 28 ** args.maximum
    if modulus_product <= dimension_bound:
        raise SystemExit("CRT modulus product does not exceed the character bound")

    moments = [
        crt([rows[index] for _, rows in parsed], moduli)
        for index in range(args.maximum + 1)
    ]
    for index, value in enumerate(moments):
        if value < 0 or value > 28 ** index:
            raise SystemExit(f"reconstructed moment m_{index} violates 0 <= m_j <= 28^j")

    first_summary = max(0, args.maximum - 12)
    expected = parse_exact_summary(args.reconstruction_log, first_summary, args.maximum)
    repeated = parse_exact_summary(
        args.reconstruction_check_log, first_summary, args.maximum
    )
    if expected != repeated:
        raise SystemExit("the two exact reconstruction summaries disagree")
    if any(moments[index] != value for index, value in expected.items()):
        raise SystemExit("CRT reconstruction disagrees with the archived exact summary")

    cross_checks = parse_source_cross_check(args.source_cross_check_log, args.maximum)
    if any(moments[index] != value for index, value in cross_checks.items()):
        raise SystemExit("CRT reconstruction disagrees with a full Pieri source moment")
    stable = read_stable(args.stable_moments, args.maximum)

    print("D4 modular Weyl exact-moment replay", flush=True)
    print("arithmetic=exact_integer_crt", flush=True)
    for path, modulus in zip(args.mod_log, moduli, strict=True):
        print(
            f"mod_log={path} sha256={digest(path)} prime={modulus} rows={args.maximum + 1}",
            flush=True,
        )
    print(
        f"reconstruction_log={args.reconstruction_log} sha256={digest(args.reconstruction_log)}",
        flush=True,
    )
    print(
        "reconstruction_check_log="
        f"{args.reconstruction_check_log} sha256={digest(args.reconstruction_check_log)}",
        flush=True,
    )
    for path in args.source_cross_check_log:
        print(f"source_cross_check_log={path} sha256={digest(path)}", flush=True)
    print(f"modulus_product={modulus_product}", flush=True)
    print(f"dimension_bound=28^{args.maximum}={dimension_bound}", flush=True)
    for index, value in enumerate(moments):
        delta = value - stable[index]
        print(
            f"D_4 exact_Weyl m_{index} = {value}; Delta_{index} = {delta}",
            flush=True,
        )
    print(
        f"SUMMARY primes={len(moduli)} moments={len(moments)} "
        f"cross_checks={len(cross_checks)} failures=0",
        flush=True,
    )
    print("__EXIT_STATUS=0", flush=True)


if __name__ == "__main__":
    main()
