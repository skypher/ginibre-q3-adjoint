#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=117.

This checks the arithmetic part of the DCT-RectTrig(E8, n=117) leaf.  The
region is a 1/20 by 1/20 rectangular union inside the root angle cap
|alpha.v| <= 4.  The local character estimates use rational cubic
minorants/majorants in w=z^2, certified against Taylor polynomials by exact
Sturm root counts on 0 <= w <= 16.

For the radial Gaussian and Weyl sine density losses, the checker uses the
fractional-part exponential lower bound

    exp(-x) >= (1457/536)^(-floor(x))
              * (1 - r + r^2/2 - r^3/6 + r^4/24
                 - r^5/120 + r^6/720 - r^7/5040),

where r = x - floor(x).  The script compares the resulting
rectangular-union lower bound for
D_E8(117) = Q_3^E8(119) - 4 Q_3^E8(117) with the negative-region bound
520*16^117.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, log10_fraction
from direct_chain_rect_e8_n131_certificate import root_identity_check


GROUP = "E8"
N = 117
D = 248
C_NEG = 8
RANK = 8
ALPHA = 250
SHAPE = 124
QUARTIC_DENOM = 50
SEXTIC_DENOM = 1800

DELTA = Fraction(4)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
S_END = Fraction(283, 5)
T_START = Fraction(30)
T_END = Fraction(50)

PI_UPPER = Fraction(355, 113)
E_UPPER = Fraction(1457, 536)
SINE_EXP_DENOM = Fraction(81, 4)

GAP_LOWER_COEFFS = (
    Fraction(99119225, 100000000),
    Fraction(-7802858, 100000000),
    Fraction(181097, 100000000),
)
GAP_UPPER_COEFFS = (
    Fraction(1),
    Fraction(-8228444, 100000000),
    Fraction(233397, 100000000),
)
COS_LOWER_COEFFS = (
    Fraction(199909270, 100000000),
    Fraction(-99730681, 100000000),
    Fraction(8092566, 100000000),
    Fraction(-216849, 100000000),
)


def trim(poly: list[Fraction]) -> list[Fraction]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_eval(poly: list[Fraction], x: Fraction) -> Fraction:
    value = Fraction(0)
    for coeff in reversed(poly):
        value = value * x + coeff
    return value


def poly_derivative(poly: list[Fraction]) -> list[Fraction]:
    return [i * poly[i] for i in range(1, len(poly))]


def poly_divmod(num: list[Fraction], den: list[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    num = trim(num[:])
    den = trim(den[:])
    if not den:
        raise ZeroDivisionError("polynomial division by zero")
    if len(num) < len(den):
        return [], num
    quotient = [Fraction(0)] * (len(num) - len(den) + 1)
    while len(num) >= len(den) and num:
        shift = len(num) - len(den)
        coeff = num[-1] / den[-1]
        quotient[shift] = coeff
        for i, den_coeff in enumerate(den):
            num[i + shift] -= coeff * den_coeff
        trim(num)
    return trim(quotient), trim(num)


def sturm_sequence(poly: list[Fraction]) -> list[list[Fraction]]:
    first = trim(poly[:])
    second = trim(poly_derivative(first))
    sequence = [first, second]
    while sequence[-1]:
        _, remainder = poly_divmod(sequence[-2], sequence[-1])
        sequence.append([-coeff for coeff in remainder])
    return sequence[:-1]


def sign(value: Fraction) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def sign_variations(values: list[Fraction]) -> int:
    nonzero = [sign(value) for value in values if value != 0]
    return sum(1 for a, b in zip(nonzero, nonzero[1:]) if a != b)


def sturm_root_count(poly: list[Fraction], left: Fraction, right: Fraction) -> int:
    sequence = sturm_sequence(poly)
    left_values = [poly_eval(item, left) for item in sequence]
    right_values = [poly_eval(item, right) for item in sequence]
    return sign_variations(left_values) - sign_variations(right_values)


def positive_no_roots(poly: list[Fraction], left: Fraction, right: Fraction) -> bool:
    return (
        poly_eval(poly, left) > 0
        and poly_eval(poly, right) > 0
        and sturm_root_count(poly, left, right) == 0
    )


def scalar_polynomial_checks() -> dict[str, bool]:
    # T12 lower for 2(1-cos sqrt(w)), divided by w.
    gap_lower_residual = [
        Fraction(0),
        Fraction(1) - GAP_LOWER_COEFFS[0],
        Fraction(-1, 12) - GAP_LOWER_COEFFS[1],
        Fraction(1, 360) - GAP_LOWER_COEFFS[2],
        Fraction(-1, 20160),
        Fraction(1, 1814400),
        Fraction(-1, 239500800),
    ]
    gap_lower_quotient = gap_lower_residual[1:]

    # Cubic upper for 2(1-cos sqrt(w)) against T10, divided by w^2.
    gap_upper_residual = [
        Fraction(0),
        GAP_UPPER_COEFFS[0] - 1,
        GAP_UPPER_COEFFS[1] + Fraction(1, 12),
        GAP_UPPER_COEFFS[2] - Fraction(1, 360),
        Fraction(1, 20160),
        Fraction(-1, 1814400),
    ]
    gap_upper_quotient = gap_upper_residual[2:]

    # T14 lower for 2 cos sqrt(w).
    cos_lower_residual = [
        Fraction(2) - COS_LOWER_COEFFS[0],
        Fraction(-1) - COS_LOWER_COEFFS[1],
        Fraction(1, 12) - COS_LOWER_COEFFS[2],
        Fraction(-1, 360) - COS_LOWER_COEFFS[3],
        Fraction(1, 20160),
        Fraction(-1, 1814400),
        Fraction(1, 239500800),
        Fraction(-1, 43589145600),
    ]

    return {
        "gap_lower_sturm": positive_no_roots(gap_lower_quotient, Fraction(0), Fraction(16)),
        "gap_upper_sturm": positive_no_roots(gap_upper_quotient, Fraction(0), Fraction(16)),
        "cos_lower_sturm": positive_no_roots(cos_lower_residual, Fraction(0), Fraction(16)),
    }


def cap_value() -> Fraction:
    # E8 has kappa = 30, so S_delta(n) = delta^2 * 30n/(4d).
    return DELTA * DELTA * 30 * N / (4 * D)


def e_upper_check() -> bool:
    # e < sum_{k=0}^7 1/k! + 1/(7!*7) < 1457/536.
    partial = sum(Fraction(1, factorial(k)) for k in range(8))
    return partial + Fraction(1, factorial(7) * 7) < E_UPPER


def exp_fractional_septic_lower(value: Fraction) -> Fraction:
    whole = value.numerator // value.denominator
    frac = value - whole
    poly = (
        Fraction(1)
        - frac
        + frac * frac / 2
        - frac**3 / 6
        + frac**4 / 24
        - frac**5 / 120
        + frac**6 / 720
        - frac**7 / 5040
    )
    return poly / (E_UPPER**whole)


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    lower = (
        GAP_LOWER_COEFFS[0] * s0
        + GAP_LOWER_COEFFS[1] * Fraction(2 * D, QUARTIC_DENOM) * s0 * s0 / N
        + GAP_LOWER_COEFFS[2] * Fraction((2 * D) ** 2, SEXTIC_DENOM) * s0**3 / (N * N)
    )
    upper = (
        GAP_UPPER_COEFFS[0] * t1
        + GAP_UPPER_COEFFS[1] * Fraction(2 * D, QUARTIC_DENOM) * t1 * t1 / N
        + GAP_UPPER_COEFFS[2] * Fraction((2 * D) ** 2, SEXTIC_DENOM) * t1**3 / (N * N)
    )
    return lower - upper


def local_h_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    return (
        (16 + 240 * COS_LOWER_COEFFS[0]) / (2 * D)
        + COS_LOWER_COEFFS[1] * (s1 + t1) / N
        + COS_LOWER_COEFFS[2]
        * Fraction(2 * D, QUARTIC_DENOM)
        * (s0 * s0 + t0 * t0)
        / (N * N)
        + COS_LOWER_COEFFS[3]
        * Fraction((2 * D) ** 2, SEXTIC_DENOM)
        * (s1**3 + t1**3)
        / (N**3)
    )


def local_direct_factor_lower(h: Fraction) -> Fraction:
    return Fraction((2 * D) ** 2) * h * h - 4


def candidate_cells() -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    cells: list[tuple[Fraction, Fraction, Fraction, Fraction]] = []
    cap = cap_value()
    s0 = S_START
    while s0 + GRID_STEP <= cap and s0 + GRID_STEP <= S_END:
        t0 = T_START
        while t0 + GRID_STEP <= T_END and t0 + GRID_STEP <= cap:
            s1 = s0 + GRID_STEP
            t1 = t0 + GRID_STEP
            if gap_lower(s0, s1, t0, t1) > 0:
                h = local_h_lower(s0, s1, t0, t1)
                if h > 0 and local_direct_factor_lower(h) > 0:
                    cells.append((s0, s1, t0, t1))
            t0 += GRID_STEP
        s0 += GRID_STEP
    return cells


def gamma_inv_square() -> Fraction:
    # Gamma(124) = 123!.
    return Fraction(1, factorial(SHAPE - 1) ** 2)


def rectangle_fraction_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction
) -> Fraction:
    s_integral = (s1**SHAPE - s0**SHAPE) / SHAPE
    t_integral = (t1**SHAPE - t0**SHAPE) / SHAPE
    return (
        gap_lower(s0, s1, t0, t1) ** 2
        * exp_fractional_septic_lower(s1 + t1)
        * s_integral
        * t_integral
        * gamma_inv_square()
        / SHAPE
    )


def a0_lower() -> Fraction:
    return a0_prefactor(ROOT_DATA[GROUP]) / (PI_UPPER**RANK)


def sine_density_lower(s1: Fraction, t1: Fraction) -> Fraction:
    exponent = Fraction(4 * D, N) * (s1 + t1) / SINE_EXP_DENOM
    return exp_fractional_septic_lower(exponent)


def negative_bound() -> int:
    return 2 * ((2 * C_NEG) ** 2 + 4) * (2 * C_NEG) ** N


def positive_rectangle_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, a0: Fraction
) -> Fraction:
    h = local_h_lower(s0, s1, t0, t1)
    direct = local_direct_factor_lower(h)
    if h <= 0 or direct <= 0:
        return Fraction(0)
    return (
        a0
        * Fraction((2 * D) ** N, N**ALPHA)
        * rectangle_fraction_lower(s0, s1, t0, t1)
        * h**N
        * direct
        * sine_density_lower(s1, t1)
    )


def main() -> None:
    cells = candidate_cells()
    if len(cells) != 111427:
        raise AssertionError(f"expected 111427 retained grid cells, got {len(cells)}")

    scalar_checks = scalar_polynomial_checks()
    gaps: list[Fraction] = []
    h_values: list[Fraction] = []
    direct_values: list[Fraction] = []
    rect_fracs: list[Fraction] = []
    positive_terms: list[Fraction] = []
    sine_exponents: list[Fraction] = []
    a0 = a0_lower()
    neg = negative_bound()

    for cell in cells:
        s0, s1, t0, t1 = cell
        gap = gap_lower(s0, s1, t0, t1)
        h = local_h_lower(s0, s1, t0, t1)
        direct = local_direct_factor_lower(h)
        rect_frac = rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = Fraction(4 * D, N) * (s1 + t1) / SINE_EXP_DENOM
        sine = exp_fractional_septic_lower(sine_exponent)
        term = (
            a0
            * Fraction((2 * D) ** N, N**ALPHA)
            * rect_frac
            * h**N
            * direct
            * sine
        )
        gaps.append(gap)
        h_values.append(h)
        direct_values.append(direct)
        rect_fracs.append(rect_frac)
        positive_terms.append(term)
        sine_exponents.append(sine_exponent)

    ratio117 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg

    checks = {
        **scalar_checks,
        "root_identity_constants": root_identity_check(),
        "e_upper_1457_536": e_upper_check(),
        "cell_count_111427": len(cells) == 111427,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h117": min(h_values) > 0,
        "positive_direct117": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio117_gt_1": ratio117 > 1,
    }

    print(f"group={GROUP} n={N} delta={DELTA} e_upper={E_UPPER}")
    print(
        f"grid_step={GRID_STEP} s_start={S_START} s_end={S_END} "
        f"t_start={T_START} t_end={T_END}"
    )
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap117_min={min(gaps)}")
    print(f"h117_min={min(h_values)} direct_factor117_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio117_lower={log10_fraction(ratio117):.2f}")
    print(f"log10_best_cell_ratio117={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=117 rectangular certificate failed")


if __name__ == "__main__":
    main()
