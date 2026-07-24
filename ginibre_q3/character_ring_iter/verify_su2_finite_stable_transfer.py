#!/usr/bin/env python3
"""Exact bounded regression for the SU(2)_k stable-transfer theorem."""

from __future__ import annotations

import argparse
from itertools import combinations_with_replacement
from typing import Dict, Iterable, Tuple

State = Dict[Tuple[int, int], int]


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-total", type=int, default=10)
    parser.add_argument("--max-factors", type=int, default=6)
    parser.add_argument("--max-label", type=int, default=8)
    args = parser.parse_args()
    if args.max_total < 0 or args.max_factors < 0 or args.max_label < 1:
        raise SystemExit("bounds must be nonnegative, with max-label at least one")

    cases = 0
    for length in range(args.max_factors + 1):
        for labels in combinations_with_replacement(
            range(1, args.max_label + 1), length
        ):
            total = sum(labels)
            if total > args.max_total:
                continue
            for minus_mask in range(1 << length):
                ordinary = signed_product(labels, minus_mask, None)
                for level in {max(1, total), max(1, total + 1)}:
                    finite = signed_product(labels, minus_mask, level)
                    if finite != ordinary:
                        raise AssertionError(
                            "stable-transfer mismatch: "
                            f"labels={labels} mask={minus_mask} level={level}"
                        )
                cases += 1

    print(
        "SU2_STABLE_TRANSFER PASS "
        f"cases={cases} max_total={args.max_total} "
        f"max_factors={args.max_factors} max_label={args.max_label} "
        "levels=T,T+1"
    )


if __name__ == "__main__":
    main()
