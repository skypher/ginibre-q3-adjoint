#!/usr/bin/env python3
"""Exact bounded partial-character regression for the SU(2)_5 fusion ring."""

from __future__ import annotations

import argparse
from itertools import product
from typing import Dict, Iterable, Tuple

State = Dict[Tuple[int, int], int]


def outputs(left: int, right: int) -> Iterable[int]:
    upper = min(left + right, 10 - left - right)
    return range(abs(left - right), upper + 1, 2)


def signed_product(exponents: Tuple[int, ...], signs: Tuple[int, ...]) -> State:
    state: State = {(0, 0): 1}
    for label in range(1, 6):
        for _ in range(exponents[label - 1]):
            next_state: State = {}
            sign = signs[label - 1]
            for (left, right), value in state.items():
                for target in outputs(left, label):
                    key = (target, right)
                    next_state[key] = next_state.get(key, 0) + value
                for target in outputs(right, label):
                    key = (left, target)
                    next_state[key] = next_state.get(key, 0) + sign * value
            state = {key: value for key, value in next_state.items() if value}
    return state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-power", type=int, default=6)
    args = parser.parse_args()
    if args.max_power < 1:
        raise SystemExit("max-power must be positive")

    cases = 0
    coefficients = 0
    minimum = None
    minimum_positive = None
    for status in product((0, 1, -1), repeat=5):
        if -1 not in status:
            continue
        active = [index for index, value in enumerate(status) if value]
        for powers in product(range(1, args.max_power + 1), repeat=len(active)):
            exponents = [0] * 5
            signs = [1] * 5
            for index, power in zip(active, powers):
                exponents[index] = power
                signs[index] = status[index]
            if sum(
                exponents[index]
                for index, value in enumerate(status)
                if value == -1
            ) % 2:
                continue
            state = signed_product(tuple(exponents), tuple(signs))
            cases += 1
            for target in range(6):
                value = state.get((target, 0), 0)
                coefficients += 1
                if value < 0:
                    raise AssertionError(
                        f"negative coefficient status={status} "
                        f"exponents={tuple(exponents)} target={target} value={value}"
                    )
                minimum = value if minimum is None else min(minimum, value)
                if value > 0:
                    minimum_positive = (
                        value if minimum_positive is None
                        else min(minimum_positive, value)
                    )

    print(
        "SU2_LEVEL_5_REGRESSION PASS "
        f"cases={cases} partial_coefficients={coefficients} "
        f"minimum={minimum} minimum_positive={minimum_positive} "
        f"max_power={args.max_power}"
    )


if __name__ == "__main__":
    main()
