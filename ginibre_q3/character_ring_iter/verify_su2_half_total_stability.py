#!/usr/bin/env python3
"""Exact bounded regression for the SU(2)_k half-total stability theorem."""

from __future__ import annotations

import argparse
from itertools import combinations_with_replacement
from typing import Dict, Iterable, Tuple

State = Dict[Tuple[int, int], int]
Mismatch = Tuple[Tuple[int, ...], int, int, int, int, int]


def ordinary_outputs(left: int, right: int) -> Iterable[int]:
    return range(abs(left - right), left + right + 1, 2)


def finite_outputs(level: int, left: int, right: int) -> Iterable[int]:
    upper = min(left + right, 2 * level - left - right)
    return range(abs(left - right), upper + 1, 2)


def signed_product(labels: Tuple[int, ...], minus_mask: int, level: int | None) -> State:
    current: State = {(0, 0): 1}
    for index, label in enumerate(labels):
        minus = bool((minus_mask >> index) & 1)
        next_state: State = {}
        for (left, right), value in current.items():
            left_outputs = (
                ordinary_outputs(left, label)
                if level is None
                else finite_outputs(level, left, label)
            )
            for output in left_outputs:
                key = (output, right)
                next_state[key] = next_state.get(key, 0) + value

            right_outputs = (
                ordinary_outputs(right, label)
                if level is None
                else finite_outputs(level, right, label)
            )
            signed_value = -value if minus else value
            for output in right_outputs:
                key = (left, output)
                next_state[key] = next_state.get(key, 0) + signed_value

        current = {key: value for key, value in next_state.items() if value}
    return current


def stable_level(labels: Tuple[int, ...], left_target: int, right_target: int = 0) -> int:
    total = sum(labels)
    largest = max(labels, default=0)
    target = max(left_target, right_target)
    return max(1, largest, left_target, right_target, (total + target + 1) // 2)


def format_mismatch(prefix: str, mismatch: Mismatch) -> str:
    labels, mask, target, level, expected, actual = mismatch
    return (
        f"{prefix} labels={labels} mask={mask} target={target} level={level} "
        f"ordinary={expected} finite={actual}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-total", type=int, default=10)
    parser.add_argument("--max-factors", type=int, default=6)
    parser.add_argument("--max-label", type=int, default=8)
    args = parser.parse_args()
    if args.max_total < 0 or args.max_factors < 0 or args.max_label < 1:
        raise SystemExit("bounds must be nonnegative, with max-label at least one")

    words = 0
    coefficient_cases = 0
    below_threshold_tests = 0
    below_threshold_mismatches = 0
    first_below_threshold_mismatch: Mismatch | None = None
    first_corner_below_threshold_mismatch: Mismatch | None = None

    for length in range(args.max_factors + 1):
        for labels in combinations_with_replacement(
            range(1, args.max_label + 1), length
        ):
            total = sum(labels)
            if total > args.max_total:
                continue
            for minus_mask in range(1 << length):
                ordinary = signed_product(labels, minus_mask, None)
                for target in range(total + 1):
                    expected = ordinary.get((target, 0), 0)
                    level = stable_level(labels, target)
                    for tested_level in {level, level + 1}:
                        finite = signed_product(labels, minus_mask, tested_level)
                        actual = finite.get((target, 0), 0)
                        if actual != expected:
                            raise AssertionError(
                                "half-total stability mismatch: "
                                f"labels={labels} mask={minus_mask} target={target} "
                                f"level={tested_level} expected={expected} actual={actual}"
                            )
                    coefficient_cases += 1

                    below = level - 1
                    largest = max(labels, default=0)
                    if below >= max(1, largest, target):
                        finite_below = signed_product(labels, minus_mask, below)
                        actual_below = finite_below.get((target, 0), 0)
                        below_threshold_tests += 1
                        if actual_below != expected:
                            mismatch = (
                                labels,
                                minus_mask,
                                target,
                                below,
                                expected,
                                actual_below,
                            )
                            below_threshold_mismatches += 1
                            if first_below_threshold_mismatch is None:
                                first_below_threshold_mismatch = mismatch
                            if target == 0 and first_corner_below_threshold_mismatch is None:
                                first_corner_below_threshold_mismatch = mismatch
                words += 1

    print(
        "SU2_HALF_TOTAL_STABLE_TRANSFER PASS "
        f"words={words} coefficient_cases={coefficient_cases} "
        f"max_total={args.max_total} max_factors={args.max_factors} "
        f"max_label={args.max_label} levels=K,K+1 "
        f"below_threshold_tests={below_threshold_tests} "
        f"below_threshold_mismatches={below_threshold_mismatches}"
    )
    if first_below_threshold_mismatch is not None:
        print(format_mismatch(
            "FIRST_BELOW_THRESHOLD_MISMATCH",
            first_below_threshold_mismatch,
        ))
    if first_corner_below_threshold_mismatch is not None:
        print(format_mismatch(
            "FIRST_CORNER_BELOW_THRESHOLD_MISMATCH",
            first_corner_below_threshold_mismatch,
        ))


if __name__ == "__main__":
    main()
