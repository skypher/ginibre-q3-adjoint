#!/usr/bin/env python3
"""Direct finite-fusion regression for SU(2) levels one through four."""

from __future__ import annotations

import argparse
from itertools import product
from typing import Dict, Iterable, Tuple

State = Dict[Tuple[int, int], int]


def outputs(level: int, left: int, right: int) -> Iterable[int]:
    upper = min(left + right, 2 * level - left - right)
    return range(abs(left - right), upper + 1, 2)


def signed_product(
    level: int,
    exponents: Tuple[int, ...],
    signs: Tuple[int, ...],
) -> State:
    state: State = {(0, 0): 1}
    for label in range(1, level + 1):
        for _ in range(exponents[label - 1]):
            next_state: State = {}
            sign = signs[label - 1]
            for (left, right), value in state.items():
                for target in outputs(level, left, label):
                    key = (target, right)
                    next_state[key] = next_state.get(key, 0) + value
                for target in outputs(level, right, label):
                    key = (left, target)
                    next_state[key] = next_state.get(key, 0) + sign * value
            state = {key: value for key, value in next_state.items() if value}
    return state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-power", type=int, default=8)
    args = parser.parse_args()
    if args.max_power < 1:
        raise SystemExit("max-power must be positive")

    grand_cases = 0
    grand_coefficients = 0
    for level in range(1, 5):
        cases = 0
        coefficients = 0
        minimum = None
        minimum_positive = None
        for status in product((0, 1, -1), repeat=level):
            if -1 not in status:
                continue
            active = [index for index, value in enumerate(status) if value]
            for powers in product(range(1, args.max_power + 1), repeat=len(active)):
                exponents = [0] * level
                signs = [1] * level
                for index, power in zip(active, powers):
                    exponents[index] = power
                    signs[index] = status[index]
                minus_parity = sum(
                    exponents[index]
                    for index, value in enumerate(status)
                    if value == -1
                )
                if minus_parity % 2:
                    continue
                state = signed_product(level, tuple(exponents), tuple(signs))
                cases += 1
                for target in range(level + 1):
                    value = state.get((target, 0), 0)
                    coefficients += 1
                    if value < 0:
                        raise AssertionError(
                            f"negative coefficient level={level} status={status} "
                            f"exponents={tuple(exponents)} target={target} value={value}"
                        )
                    minimum = value if minimum is None else min(minimum, value)
                    if value > 0:
                        minimum_positive = (
                            value if minimum_positive is None
                            else min(minimum_positive, value)
                        )
        print(
            f"SU2_LEVEL_{level}_REGRESSION PASS cases={cases} "
            f"partial_coefficients={coefficients} minimum={minimum} "
            f"minimum_positive={minimum_positive} max_power={args.max_power}"
        )
        grand_cases += cases
        grand_coefficients += coefficients
    print(
        "SU2_LEVELS_1_4_REGRESSION PASS "
        f"cases={grand_cases} partial_coefficients={grand_coefficients} "
        f"max_power={args.max_power}"
    )


if __name__ == "__main__":
    main()
