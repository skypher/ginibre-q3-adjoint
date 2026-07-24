#!/usr/bin/env python3
"""Exact bounded scan of the first non-Hall rank-six orbit chamber."""

from argparse import ArgumentParser


def outputs(left: int, right: int) -> range:
    return range(abs(left - right), min(left + right, 11 - left - right) + 1)


def step(
    state: dict[tuple[int, int], int],
    label: int,
    sign: int,
) -> dict[tuple[int, int], int]:
    result: dict[tuple[int, int], int] = {}
    for (left, right), coefficient in state.items():
        for output in outputs(left, label):
            key = (output, right)
            result[key] = result.get(key, 0) + coefficient
        for output in outputs(right, label):
            key = (left, output)
            result[key] = result.get(key, 0) + sign * coefficient
    return {key: value for key, value in result.items() if value}


def paired_step(
    state: dict[tuple[int, int], int],
    label: int,
    sign: int,
) -> dict[tuple[int, int], int]:
    return step(step(state, label, sign), label, sign)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--max-p", type=int, default=100)
    parser.add_argument("--max-q", type=int, default=100)
    args = parser.parse_args()
    if args.max_p < 0 or args.max_q < 0:
        raise SystemExit("bounds must be nonnegative")

    # Chamber: (B_1^- )^(2+2p) (B_5^+ )^(1+2q).
    base: dict[tuple[int, int], int] = {(0, 0): 1}
    base = paired_step(base, 1, -1)
    base = step(base, 5, 1)

    checks = 0
    minimum: int | None = None
    zeros: list[tuple[int, int]] = []
    p_state = base
    for p in range(args.max_p + 1):
        q_state = p_state
        for q in range(args.max_q + 1):
            value = q_state.get((0, 0), 0)
            checks += 1
            if value < 0:
                raise AssertionError(f"negative value at p={p}, q={q}: {value}")
            minimum = value if minimum is None else min(minimum, value)
            if value == 0:
                zeros.append((p, q))
            q_state = paired_step(q_state, 5, 1)
        p_state = paired_step(p_state, 1, -1)

    print(
        "SU2_O11_FIRST_CHAMBER PASS "
        f"checks={checks} max_p={args.max_p} max_q={args.max_q} "
        f"minimum={minimum} zeros={zeros}"
    )


if __name__ == "__main__":
    main()
