#!/usr/bin/env python3
"""Exact formal certificate for the SU(3), SU(4), and SU(5) recurrences.

The checker uses only Python integers.  It verifies the Bessel-determinant
reductions, carries the claimed adjoint recurrences through the binomial
transform in the normal-ordered Weyl algebra Q<w,D>/(Dw-wD-1), and reduces
the resulting differential operators modulo t*y''+y'-y=0.
"""

from __future__ import annotations

from itertools import permutations
from math import comb, factorial


# A Laurent monomial is t^(t2/2) y^yp z^zp.  Half-powers are needed only
# while checking the Bessel determinant; all differentiated expressions have
# even t2.  Coefficients remain exact integers throughout.
Monomial = tuple[int, int, int]
Laurent = dict[Monomial, int]
KPolynomial = list[int]
Weyl = dict[tuple[int, int], int]  # z^a D^b, in normal order.


def clean(value: Laurent) -> Laurent:
    return {monomial: coefficient for monomial, coefficient in value.items() if coefficient}


def add(left: Laurent, right: Laurent, right_scale: int = 1) -> Laurent:
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + right_scale * coefficient
    return clean(out)


def scale(value: Laurent, coefficient: int) -> Laurent:
    return clean({monomial: coefficient * entry for monomial, entry in value.items()})


def multiply(left: Laurent, right: Laurent) -> Laurent:
    out: Laurent = {}
    for (ta, ya, za), ca in left.items():
        for (tb, yb, zb), cb in right.items():
            monomial = (ta + tb, ya + yb, za + zb)
            out[monomial] = out.get(monomial, 0) + ca * cb
    return clean(out)


def t_shift(value: Laurent, twice_exponent: int) -> Laurent:
    return {(te + twice_exponent, ye, ze): coefficient for (te, ye, ze), coefficient in value.items()}


def determinant(matrix: list[list[Laurent]]) -> Laurent:
    size = len(matrix)
    out: Laurent = {}
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(size)
            for j in range(i + 1, size)
        )
        term: Laurent = {(0, 0, 0): -1 if inversions & 1 else 1}
        for row, column in enumerate(permutation):
            term = multiply(term, matrix[row][column])
        out = add(out, term)
    return out


def bessel_determinant(rank: int) -> Laurent:
    y: Laurent = {(0, 1, 0): 1}
    root_t_z: Laurent = {(1, 0, 1): 1}
    bessel = [y, root_t_z]
    # I_{r+1}(2 sqrt(t)) = I_{r-1}(2 sqrt(t)) - r I_r(2 sqrt(t))/sqrt(t).
    for r in range(1, rank - 1):
        bessel.append(add(bessel[r - 1], t_shift(scale(bessel[r], r), -1), -1))
    matrix = [[bessel[abs(row - column)] for column in range(rank)] for row in range(rank)]
    return determinant(matrix)


def expected_determinant(rank: int) -> Laurent:
    y: Laurent = {(0, 1, 0): 1}
    z: Laurent = {(0, 0, 1): 1}
    if rank == 3:
        # z(2y^2-yz-2tz^2)
        return {
            (0, 2, 1): 2,
            (0, 1, 2): -1,
            (2, 0, 3): -2,
        }
    if rank == 4:
        # -(4t^2z^4-8ty^2z^2+8tyz^3-tz^4+4y^4-8y^3z+4y^2z^2)/t
        return {
            (2, 0, 4): -4,
            (0, 2, 2): 8,
            (0, 1, 3): -8,
            (0, 0, 4): 1,
            (-2, 4, 0): -4,
            (-2, 3, 1): 8,
            (-2, 2, 2): -4,
        }
    if rank == 5:
        first = add(
            add(scale(multiply(y, y), 2), scale(multiply(y, z), -1)),
            {(2, 0, 2): -2, (0, 0, 2): -1},
        )
        second = {
            (2, 1, 2): -4,
            (2, 0, 3): 3,
            (0, 3, 0): 4,
            (0, 2, 1): -5,
            (0, 1, 2): 1,
        }
        return t_shift(scale(multiply(first, second), -4), -4)
    raise ValueError(rank)


def derivative(value: Laurent) -> Laurent:
    out: Laurent = {}
    for (t2, yp, zp), coefficient in value.items():
        if t2 & 1:
            raise RuntimeError("attempted to differentiate an uncancelled half-power")
        t_power = t2 // 2
        if t_power:
            monomial = (t2 - 2, yp, zp)
            out[monomial] = out.get(monomial, 0) + coefficient * t_power
        if yp:
            monomial = (t2, yp - 1, zp + 1)
            out[monomial] = out.get(monomial, 0) + coefficient * yp
        if zp:
            # z'=(y-z)/t, the reduction of t*y''+y'-y=0.
            monomial_y = (t2 - 2, yp + 1, zp - 1)
            monomial_z = (t2 - 2, yp, zp)
            out[monomial_y] = out.get(monomial_y, 0) + coefficient * zp
            out[monomial_z] = out.get(monomial_z, 0) - coefficient * zp
    return clean(out)


def theta(value: Laurent) -> Laurent:
    return t_shift(derivative(value), 2)


def sequence_shift(value: Laurent) -> Laurent:
    # S=(theta+1)D sends sum a_k t^k/(k!)^2 to
    # sum a_{k+1} t^k/(k!)^2.
    first = derivative(value)
    return add(theta(first), first)


def kpoly_clean(value: KPolynomial) -> KPolynomial:
    while len(value) > 1 and value[-1] == 0:
        value.pop()
    return value


def kpoly_add(left: KPolynomial, right: KPolynomial, right_scale: int = 1) -> KPolynomial:
    out = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        out[index] += coefficient
    for index, coefficient in enumerate(right):
        out[index] += right_scale * coefficient
    return kpoly_clean(out)


def kpoly_multiply(left: KPolynomial, right: KPolynomial) -> KPolynomial:
    out = [0] * (len(left) + len(right) - 1)
    for i, ca in enumerate(left):
        for j, cb in enumerate(right):
            out[i + j] += ca * cb
    return kpoly_clean(out)


def kpoly_evaluate(value: KPolynomial, argument: int) -> int:
    out = 0
    for coefficient in reversed(value):
        out = out * argument + coefficient
    return out


def apply_theta_polynomial(polynomial: KPolynomial, value: Laurent) -> Laurent:
    out: Laurent = {}
    for coefficient in reversed(polynomial):
        out = add(theta(out), scale(value, coefficient))
    return out


def weyl_add(left: Weyl, right: Weyl) -> Weyl:
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def falling(number: int, length: int) -> int:
    out = 1
    for offset in range(length):
        out *= number - offset
    return out


def weyl_multiply(left: Weyl, right: Weyl) -> Weyl:
    out: Weyl = {}
    # D^b z^c = sum_l binom(b,l) (c)_l z^(c-l) D^(b-l).
    for (a, b), ca in left.items():
        for (c, d), cb in right.items():
            for ell in range(min(b, c) + 1):
                monomial = (a + c - ell, b + d - ell)
                coefficient = ca * cb * comb(b, ell) * falling(c, ell)
                out[monomial] = out.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def weyl_power(value: Weyl, exponent: int) -> Weyl:
    out: Weyl = {(0, 0): 1}
    for _ in range(exponent):
        out = weyl_multiply(out, value)
    return out


def polynomial_of_weyl(polynomial: KPolynomial, value: Weyl) -> Weyl:
    out: Weyl = {}
    power: Weyl = {(0, 0): 1}
    for coefficient in polynomial:
        if coefficient:
            out = weyl_add(out, {monomial: coefficient * entry for monomial, entry in power.items()})
        power = weyl_multiply(power, value)
    return out


def transformed_operator(recurrence: list[KPolynomial]) -> tuple[int, dict[int, KPolynomial], Weyl]:
    # If M is the EGF of m and A=e^z M is the EGF of its inverse binomial
    # transform E, then D -> D-1 and theta=zD -> zD-z under conjugation.
    theta_minus_z: Weyl = {(1, 1): 1, (1, 0): -1}
    d_minus_one: Weyl = {(0, 1): 1, (0, 0): -1}
    operator: Weyl = {}
    for shift, polynomial in enumerate(recurrence):
        term = weyl_multiply(
            polynomial_of_weyl(polynomial, theta_minus_z),
            weyl_power(d_minus_one, shift),
        )
        operator = weyl_add(operator, term)

    index_shift = max(a - b for a, b in operator)
    q: dict[int, KPolynomial] = {}
    for (a, b), coefficient in operator.items():
        shift = index_shift + b - a
        falling_polynomial: KPolynomial = [1]
        for offset in range(a):
            falling_polynomial = kpoly_multiply(
                falling_polynomial, [index_shift - offset, 1]
            )
        q[shift] = kpoly_add(q.get(shift, [0]), scale_kpoly(falling_polynomial, coefficient))
    q = {shift: kpoly_clean(polynomial) for shift, polynomial in q.items() if any(polynomial)}
    return index_shift, q, operator


def scale_kpoly(value: KPolynomial, coefficient: int) -> KPolynomial:
    return [coefficient * entry for entry in value]


def partition_dimensions(total: int, maximum_length: int) -> int:
    answer = 0

    def visit(remaining: int, maximum_part: int, partition: list[int]) -> None:
        nonlocal answer
        if remaining == 0:
            length = len(partition)
            numerator = factorial(total)
            denominator = 1
            for i, part in enumerate(partition):
                denominator *= factorial(part + length - i - 1)
            for i in range(length):
                for j in range(i + 1, length):
                    numerator *= partition[i] - partition[j] + j - i
            if numerator % denominator:
                raise RuntimeError("nonintegral hook dimension")
            dimension = numerator // denominator
            answer += dimension * dimension
            return
        if len(partition) == maximum_length:
            return
        for part in range(min(maximum_part, remaining), 0, -1):
            partition.append(part)
            visit(remaining - part, part, partition)
            partition.pop()

    visit(total, total, [])
    return answer


def trace_and_adjoint_moments(rank: int, maximum: int) -> tuple[list[int], list[int]]:
    trace = [partition_dimensions(degree, rank) for degree in range(maximum + 1)]
    adjoint = []
    for degree in range(maximum + 1):
        adjoint.append(
            sum(
                (-1 if (degree - j) & 1 else 1) * comb(degree, j) * trace[j]
                for j in range(degree + 1)
            )
        )
    return trace, adjoint


def check_low_conjugated_coefficients(
    operator: Weyl, index_shift: int, trace: list[int]
) -> None:
    for degree in range(index_shift):
        coefficient = 0
        for (a, b), entry in operator.items():
            if degree < a:
                continue
            index = degree - a + b
            coefficient += entry * falling(degree, a) * trace[index]
        if coefficient:
            raise RuntimeError(
                f"nonzero omitted conjugated coefficient degree={degree}: {coefficient}"
            )


def expected_q() -> dict[int, tuple[int, dict[int, KPolynomial]]]:
    return {
        3: (
            1,
            {
                0: [9, 18, 9],
                1: [-77, -78, -19],
                2: [109, 70, 11],
                3: [-25, -10, -1],
            },
        ),
        4: (
            1,
            {
                0: [-128, -320, -256, -64],
                1: [2732, 3582, 1526, 212],
                2: [-10482, -9223, -2662, -253],
                3: [12884, 8378, 1795, 127],
                4: [-5312, -2625, -428, -23],
                5: [576, 208, 25, 1],
            },
        ),
        5: (
            2,
            {
                0: [900, 2700, 2925, 1350, 225],
                1: [-41268, -66172, -38957, -9962, -934],
                2: [297420, 324256, 130422, 22936, 1487],
                3: [-689052, -562992, -170218, -22580, -1108],
                4: [601200, 384340, 91097, 9494, 367],
                5: [-176688, -87412, -15993, -1282, -38],
                6: [14400, 5280, 724, 44, 1],
            },
        ),
    }


def recurrence_data() -> dict[int, list[KPolynomial]]:
    # For N=3 the variable is h=k-1, so the shifts are m_h,m_{h+1},m_{h+2}.
    return {
        3: [
            [16, 24, 8],
            [18, 25, 7],
            [-16, -8, -1],
        ],
        4: [
            [270, 495, 270, 45],
            [972, 1242, 522, 72],
            [-750, -340, 0, 10],
            [-1524, -1054, -230, -16],
            [392, 161, 22, 1],
        ],
        5: [
            [6912, 13824, 9024, 2304, 192],
            [29664, 42384, 21776, 4704, 352],
            [-3336, 6598, 6305, 1632, 129],
            [-30588, -23657, -6413, -732, -30],
            [6400, 2880, 484, 36, 1],
        ],
    }


def verify_rank(rank: int, recurrence: list[KPolynomial]) -> None:
    determinant_value = bessel_determinant(rank)
    expected = expected_determinant(rank)
    if determinant_value != expected:
        raise RuntimeError(f"SU({rank}) Bessel determinant reduction failed")

    index_shift, q, conjugated = transformed_operator(recurrence)
    expected_shift, expected_polynomials = expected_q()[rank]
    if index_shift != expected_shift or q != expected_polynomials:
        raise RuntimeError(
            f"SU({rank}) binomial/Weyl transform failed: shift={index_shift}, q={q}"
        )

    shifted = determinant_value
    total: Laurent = {}
    for shift in range(max(q) + 1):
        if shift in q:
            total = add(total, apply_theta_polynomial(q[shift], shifted))
        shifted = sequence_shift(shifted)
    if total:
        raise RuntimeError(f"SU({rank}) reduced differential identity is nonzero: {total}")

    trace, adjoint = trace_and_adjoint_moments(rank, 24)
    check_low_conjugated_coefficients(conjugated, index_shift, trace)
    for k in range(0, 25 - len(recurrence)):
        residual = sum(
            kpoly_evaluate(polynomial, k) * adjoint[k + shift]
            for shift, polynomial in enumerate(recurrence)
        )
        if residual:
            raise RuntimeError(f"SU({rank}) recurrence residual at k={k}: {residual}")

    print(
        f"SU({rank}) determinant_terms={len(determinant_value)} "
        f"weyl_terms={len(conjugated)} binomial_index_shift={index_shift} "
        f"q_shifts={len(q)} reduced_identity_terms=0 finite_residuals=0"
    )


def main() -> None:
    for rank, recurrence in recurrence_data().items():
        verify_rank(rank, recurrence)
    print(
        "SUMMARY SU(3..5) recurrence derivations verified by exact integer "
        "Bessel/Weyl/differential algebra: ALL PASS"
    )


if __name__ == "__main__":
    main()
