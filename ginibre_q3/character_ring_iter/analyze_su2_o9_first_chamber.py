#!/usr/bin/env python3
"""Exact integer analyzer for the first unresolved rank-five O_9 chamber.

The chamber is the post-separated-pole transport obstruction

    B_1^- exponent 2+2p,
    B_2^+ exponent 1,
    B_4^- exponent 2+2q,

with p,q >= 0.  The program uses only integral fusion matrices.  It is a
bounded diagnostic and recurrence-discovery tool, not an all-index proof.
"""

from __future__ import annotations

import argparse
from typing import Iterable, List, Sequence

Vector = List[int]
Matrix = List[List[int]]


def zero_matrix(size: int) -> Matrix:
    return [[0 for _ in range(size)] for _ in range(size)]


def identity_matrix(size: int) -> Matrix:
    result = zero_matrix(size)
    for index in range(size):
        result[index][index] = 1
    return result


def matrix_add(left: Matrix, right: Matrix, sign: int = 1) -> Matrix:
    return [
        [a + sign * b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    result = zero_matrix(size)
    for row in range(size):
        for middle, coefficient in enumerate(left[row]):
            if coefficient == 0:
                continue
            for column, value in enumerate(right[middle]):
                if value:
                    result[row][column] += coefficient * value
    return result


def matrix_vector(matrix: Matrix, vector: Sequence[int]) -> Vector:
    return [
        sum(coefficient * value for coefficient, value in zip(row, vector))
        for row in matrix
    ]


def kronecker(left: Matrix, right: Matrix) -> Matrix:
    left_size = len(left)
    right_size = len(right)
    size = left_size * right_size
    result = zero_matrix(size)
    for left_row in range(left_size):
        for left_column in range(left_size):
            left_value = left[left_row][left_column]
            if left_value == 0:
                continue
            for right_row in range(right_size):
                for right_column in range(right_size):
                    right_value = right[right_row][right_column]
                    if right_value:
                        result[left_row * right_size + right_row][
                            left_column * right_size + right_column
                        ] += left_value * right_value
    return result


def orbit_fusion_matrix(rank: int, label: int) -> Matrix:
    """Multiplication by B_label in O_(2*rank-1), columns are sources."""
    result = zero_matrix(rank)
    for source in range(rank):
        lower = abs(label - source)
        upper = min(label + source, 2 * rank - 1 - label - source)
        for target in range(lower, upper + 1):
            result[target][source] += 1
    return result


def signed_tensor_operator(fusion: Matrix, sign: int) -> Matrix:
    rank = len(fusion)
    identity = identity_matrix(rank)
    return matrix_add(
        kronecker(fusion, identity),
        kronecker(identity, fusion),
        sign=sign,
    )


def square(matrix: Matrix) -> Matrix:
    return matrix_multiply(matrix, matrix)


def apply_product(operators: Iterable[Matrix], vector: Vector) -> Vector:
    result = vector
    for operator in reversed(list(operators)):
        result = matrix_vector(operator, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-p", type=int, default=100)
    parser.add_argument("--max-q", type=int, default=100)
    parser.add_argument("--print-grid", type=int, default=8)
    args = parser.parse_args()
    if args.max_p < 0 or args.max_q < 0 or args.print_grid < 0:
        raise SystemExit("bounds must be nonnegative")

    rank = 5
    fusion = [orbit_fusion_matrix(rank, label) for label in range(rank)]
    minus_1 = signed_tensor_operator(fusion[1], -1)
    plus_2 = signed_tensor_operator(fusion[2], +1)
    minus_4 = signed_tensor_operator(fusion[4], -1)

    step_p = square(minus_1)
    step_q = square(minus_4)
    if matrix_multiply(step_p, step_q) != matrix_multiply(step_q, step_p):
        raise AssertionError("the two residual transfer operators do not commute")

    vacuum = [0] * (rank * rank)
    vacuum[0] = 1
    base = apply_product(
        [square(minus_1), plus_2, square(minus_4)],
        vacuum,
    )

    values: List[List[int]] = []
    q_vector = base
    minimum: int | None = None
    minimum_at = (0, 0)
    first_negative: tuple[int, int, int] | None = None

    for q in range(args.max_q + 1):
        row: List[int] = []
        p_vector = q_vector
        for p in range(args.max_p + 1):
            value = p_vector[0]
            row.append(value)
            if minimum is None or value < minimum:
                minimum = value
                minimum_at = (p, q)
            if value < 0 and first_negative is None:
                first_negative = (p, q, value)
            p_vector = matrix_vector(step_p, p_vector)
        values.append(row)
        q_vector = matrix_vector(step_q, q_vector)

    print(
        "SU2_O9_FIRST_CHAMBER"
        f" max_p={args.max_p} max_q={args.max_q}"
        f" minimum={minimum}@p{minimum_at[0]},q{minimum_at[1]}"
        f" result={'FAIL' if first_negative else 'PASS'}"
    )
    if first_negative is not None:
        p, q, value = first_negative
        print(f"FIRST_NEGATIVE p={p} q={q} value={value}")

    displayed_p = min(args.print_grid, args.max_p + 1)
    displayed_q = min(args.print_grid, args.max_q + 1)
    if displayed_p and displayed_q:
        print("GRID rows=q columns=p")
        for q in range(displayed_q):
            print(" ".join(str(values[q][p]) for p in range(displayed_p)))


if __name__ == "__main__":
    main()
