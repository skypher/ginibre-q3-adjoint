#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=115.

This checks the arithmetic part of the DCT-RectTrig(E8, n=115) leaf.  The
region is a 1/20 by 1/20 rectangular union inside the root angle cap
|alpha.v| <= 4.  It keeps the cubic gap and average-character scalar
polynomials from the n=117 bridge step, and replaces the linear Weyl-sine
loss by the two-term bound

    log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100.

The sine scalar inequality is certified by an exact Sturm root count on
0 <= (z/2)^2 <= 4.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    COS_LOWER_COEFFS,
    GAP_LOWER_COEFFS,
    GAP_UPPER_COEFFS,
    positive_no_roots,
    scalar_polynomial_checks,
)
from direct_chain_rect_e8_n131_certificate import root_identity_check


GROUP = "E8"
N = 115
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
S_END = Fraction(278, 5)
T_START = Fraction(30)
T_END = Fraction(48)

PI_UPPER = Fraction(355, 113)
E_UPPER = Fraction(1457, 536)
SINE_QUADRATIC = Fraction(1, 24)
SINE_QUARTIC = Fraction(1, 2100)


def poly_trim(poly: list[Fraction]) -> list[Fraction]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    result = [Fraction(0)] * size
    for i in range(size):
        result[i] = (left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0)
    return poly_trim(result)


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_coeff in enumerate(left):
        for j, right_coeff in enumerate(right):
            result[i + j] += left_coeff * right_coeff
    return poly_trim(result)


def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
    return poly_trim([factor * coeff for coeff in poly])


def poly_pow(poly: list[Fraction], exponent: int) -> list[Fraction]:
    result = [Fraction(1)]
    for _ in range(exponent):
        result = poly_mul(result, poly)
    return result


def sine_scalar_check() -> bool:
    # Let y=(z/2)^2.  The target exponent is y/6 + 4y^2/525.
    exponent_poly = [Fraction(0), Fraction(1, 6), Fraction(4, 525)]
    sinc_lower = [Fraction((-1) ** k, factorial(2 * k + 1)) for k in range(8)]
    exp_poly: list[Fraction] = []
    for k in range(6):
        exp_poly = poly_add(exp_poly, poly_scale(poly_pow(exponent_poly, k), Fraction(1, factorial(k))))
    residual = poly_add(poly_mul(sinc_lower, exp_poly), [Fraction(-1)])
    quotient = residual[:]
    while quotient and quotient[0] == 0:
        quotient = quotient[1:]
    return positive_no_roots(quotient, Fraction(0), Fraction(4))


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
    exponent = (
        2 * SINE_QUADRATIC * Fraction(2 * D) * (s1 + t1) / N
        + 2
        * SINE_QUARTIC
        * Fraction((2 * D) ** 2, QUARTIC_DENOM)
        * (s1 * s1 + t1 * t1)
        / (N * N)
    )
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
    if len(cells) != 98455:
        raise AssertionError(f"expected 98455 retained grid cells, got {len(cells)}")

    cubic_checks = scalar_polynomial_checks()
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
        sine_exponent = (
            2 * SINE_QUADRATIC * Fraction(2 * D) * (s1 + t1) / N
            + 2
            * SINE_QUARTIC
            * Fraction((2 * D) ** 2, QUARTIC_DENOM)
            * (s1 * s1 + t1 * t1)
            / (N * N)
        )
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

    ratio115 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg

    checks = {
        **cubic_checks,
        "sine_quartic_sturm": sine_scalar_check(),
        "root_identity_constants": root_identity_check(),
        "e_upper_1457_536": e_upper_check(),
        "cell_count_98455": len(cells) == 98455,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h115": min(h_values) > 0,
        "positive_direct115": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio115_gt_1": ratio115 > 1,
    }

    print(f"group={GROUP} n={N} delta={DELTA} e_upper={E_UPPER}")
    print(
        f"grid_step={GRID_STEP} s_start={S_START} s_end={S_END} "
        f"t_start={T_START} t_end={T_END}"
    )
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap115_min={min(gaps)}")
    print(f"h115_min={min(h_values)} direct_factor115_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio115_lower={log10_fraction(ratio115):.2f}")
    print(f"log10_best_cell_ratio115={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=115 rectangular certificate failed")


if __name__ == "__main__":
    main()
