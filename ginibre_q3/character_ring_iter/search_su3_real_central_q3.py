#!/usr/bin/env python3
"""Exact bounded Q3 search in real central SU(3) character generators."""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement


Exponent = tuple[int, int]
Laurent = dict[Exponent, int]


def multiply(left: Laurent, right: Laurent) -> Laurent:
    answer: Laurent = {}
    for (a, b), left_coefficient in left.items():
        for (c, d), right_coefficient in right.items():
            exponent = (a + c, b + d)
            answer[exponent] = (
                answer.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return answer


def add(left: Laurent, right: Laurent) -> Laurent:
    answer = dict(left)
    for exponent, coefficient in right.items():
        answer[exponent] = answer.get(exponent, 0) + coefficient
    return answer


@lru_cache(maxsize=None)
def irreducible_character(p: int, q: int) -> Laurent:
    """Return the SU(3) character of highest weight (p,q).

    The Gelfand--Tsetlin patterns for the GL(3) partition (p+q,q,0)
    enumerate its weight monomials.  We impose x1*x2*x3=1 and record the
    exponents of x1 and x2.
    """

    highest = (p + q, q, 0)
    answer: Laurent = {}
    for middle_left in range(highest[1], highest[0] + 1):
        for middle_right in range(highest[2], highest[1] + 1):
            for bottom in range(middle_right, middle_left + 1):
                weight_1 = bottom
                weight_2 = middle_left + middle_right - bottom
                weight_3 = sum(highest) - middle_left - middle_right
                exponent = (weight_1 - weight_3, weight_2 - weight_3)
                answer[exponent] = answer.get(exponent, 0) + 1
    return answer


def weyl_density() -> Laurent:
    answer: Laurent = {(0, 0): 1}
    roots = ((1, -1), (-1, 1), (2, 1), (-2, -1), (1, 2), (-1, -2))
    for root in roots:
        answer = multiply(answer, {(0, 0): 1, root: -1})
    return answer


WEYL_DENSITY = weyl_density()


def invariant_multiplicity(characters: tuple[Laurent, ...]) -> int:
    product: Laurent = {(0, 0): 1}
    for character in characters:
        product = multiply(product, character)
    numerator = sum(
        coefficient * product.get((-a, -b), 0)
        for (a, b), coefficient in WEYL_DENSITY.items()
    )
    if numerator % 6 != 0:
        raise AssertionError(f"nonintegral Weyl constant term {numerator}/6")
    return numerator // 6


def make_invariant_function(generators: tuple[Laurent, ...]):
    @lru_cache(maxsize=None)
    def invariant(labels: tuple[int, ...]) -> int:
        return invariant_multiplicity(tuple(generators[index] for index in labels))

    return invariant


def q3_value(
    minus: tuple[int, ...],
    plus: tuple[int, ...],
    invariant,
) -> int:
    labels = minus + plus
    answer = 0
    for mask in range(1 << len(labels)):
        right = tuple(
            sorted(labels[index] for index in range(len(labels)) if mask >> index & 1)
        )
        left = tuple(
            sorted(labels[index] for index in range(len(labels)) if not mask >> index & 1)
        )
        parity = sum(mask >> index & 1 for index in range(len(minus))) & 1
        term = invariant(left) * invariant(right)
        answer += -term if parity else term
    return answer


def search(name: str, labels: tuple[str, ...], generators: tuple[Laurent, ...]) -> int:
    invariant = make_invariant_function(generators)
    total_cases = 0
    overall_minimum: int | None = None
    for factors in range(2, 9):
        degree_cases = 0
        degree_minimum: int | None = None
        degree_witness: tuple[tuple[int, ...], tuple[int, ...]] | None = None
        for minus_count in range(2, factors + 1, 2):
            for minus in combinations_with_replacement(range(len(labels)), minus_count):
                for plus in combinations_with_replacement(
                    range(len(labels)), factors - minus_count
                ):
                    value = q3_value(minus, plus, invariant)
                    degree_cases += 1
                    if degree_minimum is None or value < degree_minimum:
                        degree_minimum = value
                        degree_witness = (minus, plus)
        if degree_minimum is None or degree_witness is None:
            raise AssertionError("empty search degree")
        total_cases += degree_cases
        if overall_minimum is None or degree_minimum < overall_minimum:
            overall_minimum = degree_minimum
        minus_names = [labels[index] for index in degree_witness[0]]
        plus_names = [labels[index] for index in degree_witness[1]]
        print(
            f"SU3_REAL_CENTRAL_Q3 suite={name} factors={factors} "
            f"cases={degree_cases} minimum={degree_minimum} "
            f"minus={minus_names} plus={plus_names}",
            flush=True,
        )
        if degree_minimum < 0:
            print("SU3_REAL_CENTRAL_Q3 SEARCH: NEGATIVE WITNESS FOUND")
            return total_cases
    print(
        f"SU3_REAL_CENTRAL_Q3 suite={name} total_cases={total_cases} "
        f"overall_minimum={overall_minimum} PASS"
    )
    return total_cases


def main() -> int:
    adjoint = irreducible_character(1, 1)
    adjoint_invariant = make_invariant_function((adjoint,))
    expected_adjoint_moments = (1, 0, 1, 2, 8, 32, 145, 702, 3598, 19280)
    actual_adjoint_moments = tuple(
        adjoint_invariant((0,) * degree)
        for degree in range(len(expected_adjoint_moments))
    )
    if actual_adjoint_moments != expected_adjoint_moments:
        raise AssertionError("SU(3) adjoint moment prefix mismatch")
    print("SU3_REAL_CENTRAL_Q3 adjoint_moment_prefix=PASS")

    self_dual_labels = tuple(f"chi_({p},{p})" for p in range(1, 5))
    self_dual_generators = tuple(irreducible_character(p, p) for p in range(1, 5))
    self_dual_cases = search(
        "self_dual_p_le_4", self_dual_labels, self_dual_generators
    )

    real_labels = (
        "chi_(1,0)+chi_(0,1)",
        "chi_(1,1)",
        "chi_(2,0)+chi_(0,2)",
        "chi_(2,1)+chi_(1,2)",
        "chi_(3,0)+chi_(0,3)",
    )
    real_generators = (
        add(irreducible_character(1, 0), irreducible_character(0, 1)),
        adjoint,
        add(irreducible_character(2, 0), irreducible_character(0, 2)),
        add(irreducible_character(2, 1), irreducible_character(1, 2)),
        add(irreducible_character(3, 0), irreducible_character(0, 3)),
    )
    real_cases = search("real_generators", real_labels, real_generators)
    print(
        "SU3_REAL_CENTRAL_Q3 SEARCH: NO NEGATIVE IN TESTED DOMAIN "
        f"total_cases={self_dual_cases + real_cases}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
