#!/usr/bin/env python3
"""Exact audit of the stable full-Q3 identities in full_q3_extension.tex.

This is a short rational replay, not a substitute for the proofs.  It checks
the two independent derivations of each exponential coefficient and compares
the stable classical moments with the vendored exact moment row.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from pathlib import Path


def moments_from_cumulants(cumulants: list[Fraction]) -> list[Fraction]:
    maximum = len(cumulants) - 1
    moments = [Fraction(0) for _ in range(maximum + 1)]
    moments[0] = Fraction(1)
    for n in range(1, maximum + 1):
        moments[n] = sum(
            Fraction(comb(n - 1, k - 1)) * cumulants[k] * moments[n - k]
            for k in range(1, n + 1)
        )
    return moments


def direct_integral(moment: list[Fraction], a: int, n: int) -> Fraction:
    return sum(
        Fraction((-1) ** j * comb(2 * a, j) * comb(n, k))
        * moment[2 * a - j + k]
        * moment[j + n - k]
        for j in range(2 * a + 1)
        for k in range(n + 1)
    )


def odd_double_factorial(index: int) -> int:
    if index == -1:
        return 1
    answer = 1
    for value in range(1, index + 1, 2):
        answer *= value
    return answer


def stable_su_moments(maximum: int) -> list[Fraction]:
    return [
        Fraction(
            sum(
                (-1) ** (degree - j) * comb(degree, j) * factorial(j)
                for j in range(degree + 1)
            )
        )
        for degree in range(maximum + 1)
    ]


def stable_su_formula(a: int, maximum_n: int) -> list[Fraction]:
    cumulants = [Fraction(0) for _ in range(maximum_n + 1)]
    if maximum_n >= 1:
        cumulants[1] = Fraction(2 * a)
    for degree in range(2, maximum_n + 1):
        cumulants[degree] = Fraction((2 * a + 2) * factorial(degree - 1))
    base = moments_from_cumulants(cumulants)
    return [Fraction(factorial(2 * a)) * value for value in base]


def stable_classical_moments(maximum: int) -> list[Fraction]:
    cumulants = [Fraction(0) for _ in range(maximum + 1)]
    if maximum >= 2:
        cumulants[2] = Fraction(1)
    for degree in range(3, maximum + 1):
        cumulants[degree] = Fraction(factorial(degree - 1), 2)
    return moments_from_cumulants(cumulants)


def h_j_moments(j: int, maximum_n: int) -> list[Fraction]:
    cumulants = [Fraction(0) for _ in range(maximum_n + 1)]
    if maximum_n >= 1:
        cumulants[1] = Fraction(2 * j)
    if maximum_n >= 2:
        cumulants[2] = Fraction(2 * j + 2)
    for degree in range(3, maximum_n + 1):
        cumulants[degree] = Fraction((2 * j + 1) * factorial(degree - 1))
    return moments_from_cumulants(cumulants)


def stable_classical_formula(a: int, maximum_n: int) -> list[Fraction]:
    answer = [Fraction(0) for _ in range(maximum_n + 1)]
    for j in range(a + 1):
        gamma = (
            Fraction(comb(2 * a, 2 * j) * odd_double_factorial(2 * a - 2 * j - 1))
            * Fraction(comb(2 * j, j), 4**j)
            * factorial(2 * j)
        )
        component = h_j_moments(j, maximum_n)
        for n in range(maximum_n + 1):
            answer[n] += gamma * component[n]
    return answer


def read_vendored_classical_moments() -> list[int]:
    path = Path(__file__).resolve().parent / "references" / "oeis_A002137_stable.txt"
    rows: list[int] = []
    expected_index = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.split()
        if not fields:
            continue
        if len(fields) != 2:
            raise RuntimeError(f"malformed stable moment row: {raw!r}")
        index, value = map(int, fields)
        if index != expected_index:
            raise RuntimeError(
                f"nonconsecutive stable moment row: got {index}, expected {expected_index}"
            )
        rows.append(value)
        expected_index += 1
    if not rows:
        raise RuntimeError("empty stable moment table")
    return rows


def main() -> None:
    maximum_a = 10
    maximum_n = 24
    maximum_degree = 2 * maximum_a + maximum_n

    su_moments = stable_su_moments(maximum_degree)
    classical_moments = stable_classical_moments(maximum_degree)

    checked = 0
    for a in range(1, maximum_a + 1):
        su_formula = stable_su_formula(a, maximum_n)
        classical_formula = stable_classical_formula(a, maximum_n)
        for n in range(maximum_n + 1):
            su_direct = direct_integral(su_moments, a, n)
            classical_direct = direct_integral(classical_moments, a, n)
            if su_direct != su_formula[n]:
                raise RuntimeError(
                    f"stable SU identity failed at a={a}, n={n}: "
                    f"direct={su_direct}, formula={su_formula[n]}"
                )
            if classical_direct != classical_formula[n]:
                raise RuntimeError(
                    f"stable classical identity failed at a={a}, n={n}: "
                    f"direct={classical_direct}, formula={classical_formula[n]}"
                )
            if su_direct <= 0 or classical_direct <= 0:
                raise RuntimeError(
                    f"nonpositive stable coefficient at a={a}, n={n}: "
                    f"SU={su_direct}, classical={classical_direct}"
                )
            checked += 2

    vendored = read_vendored_classical_moments()
    regenerated = stable_classical_moments(len(vendored) - 1)
    for index, (expected, actual) in enumerate(zip(vendored, regenerated, strict=True)):
        if actual.denominator != 1 or actual.numerator != expected:
            raise RuntimeError(
                f"vendored stable moment mismatch at degree {index}: "
                f"table={expected}, regenerated={actual}"
            )

    print("full-Q3 stable identity audit: PASS")
    print(f"exact hierarchy identities checked: {checked}")
    print(f"vendored stable moments checked: {len(vendored)}")
    print(f"parameter box: 1 <= a <= {maximum_a}, 0 <= n <= {maximum_n}")


if __name__ == "__main__":
    main()
