#!/usr/bin/env python3
"""Exact regression for odd-level SU(2) simple-current orbit lifting."""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from typing import Dict, Iterable, Tuple

Labels = Tuple[int, ...]


def outputs(level: int, left: int, right: int) -> Iterable[int]:
    upper = min(left + right, 2 * level - left - right)
    return range(abs(left - right), upper + 1, 2)


@lru_cache(maxsize=None)
def distribution(level: int, labels: Labels) -> Tuple[Tuple[int, int], ...]:
    current: Dict[int, int] = {0: 1}
    for label in labels:
        following: Dict[int, int] = defaultdict(int)
        for incoming, multiplicity in current.items():
            for outgoing in outputs(level, incoming, label):
                following[outgoing] += multiplicity
        current = following
    return tuple(sorted(current.items()))


def coefficient(level: int, labels: Labels, target: int) -> int:
    return dict(distribution(level, tuple(sorted(labels)))).get(target, 0)


def signed_corner(level: int, labels: Labels, minus_mask: int) -> int:
    length = len(labels)
    total = 0
    for subset in range(1 << length):
        left = tuple(labels[i] for i in range(length) if not ((subset >> i) & 1))
        right = tuple(labels[i] for i in range(length) if (subset >> i) & 1)
        sign = -1 if (subset & minus_mask).bit_count() & 1 else 1
        total += sign * coefficient(level, left, 0) * coefficient(level, right, 0)
    return total


def even_lift(level: int, label: int) -> int:
    return label if label % 2 == 0 else level - label


def transformed_mask(labels: Labels, minus_mask: int) -> int:
    result = minus_mask
    for index, label in enumerate(labels):
        if label % 2:
            result ^= 1 << index
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-level", type=int, default=9)
    parser.add_argument("--max-total", type=int, default=12)
    parser.add_argument("--max-factors", type=int, default=5)
    args = parser.parse_args()
    if args.max_level < 1 or args.max_level % 2 == 0:
        raise SystemExit("max-level must be a positive odd integer")

    cases = 0
    forced_zero = 0
    two_term = 0
    for level in range(1, args.max_level + 1, 2):
        level_cases = 0
        for length in range(args.max_factors + 1):
            for labels in combinations_with_replacement(range(level + 1), length):
                if sum(labels) > args.max_total:
                    continue
                lifts = tuple(even_lift(level, label) for label in labels)
                odd_count = sum(label % 2 for label in labels)
                for minus_mask in range(1 << length):
                    actual = signed_corner(level, labels, minus_mask)
                    if odd_count % 2:
                        expected = 0
                        forced_zero += 1
                    else:
                        companion_mask = transformed_mask(labels, minus_mask)
                        first = signed_corner(level, lifts, minus_mask)
                        second = signed_corner(level, lifts, companion_mask)
                        if (first + second) % 2:
                            raise AssertionError(
                                f"nonintegral average level={level} labels={labels} "
                                f"mask={minus_mask} first={first} second={second}"
                            )
                        expected = (first + second) // 2
                        two_term += 1
                    if actual != expected:
                        raise AssertionError(
                            f"lift mismatch level={level} labels={labels} "
                            f"mask={minus_mask} actual={actual} expected={expected}"
                        )
                    cases += 1
                    level_cases += 1
        print(f"ODD_SIMPLE_CURRENT_LIFT level={level} cases={level_cases} result=PASS")

    print(
        "SU2_ODD_SIMPLE_CURRENT_LIFT PASS "
        f"cases={cases} forced_zero={forced_zero} two_term={two_term} "
        f"max_level={args.max_level} max_total={args.max_total} "
        f"max_factors={args.max_factors}"
    )


if __name__ == "__main__":
    main()
