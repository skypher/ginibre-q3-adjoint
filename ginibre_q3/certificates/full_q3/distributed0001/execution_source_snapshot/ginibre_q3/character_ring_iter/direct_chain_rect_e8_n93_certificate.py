#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=93.

This replays the n=95 rectangular machinery with a slightly wider root-angle
cap |alpha.v| <= 21/5.  The local scalar bounds are rechecked on
0 <= z^2 <= (21/5)^2.  The gap lower cubic uses the adjusted coefficient
7/4000 on z^6, and the Weyl-sine estimate uses the quartic bound

    log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2000.

The negative-region estimate uses the same shifted Chebyshev degree-27
majorant as the n=95 certificate.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n99_n97_certificate as shared
import direct_chain_rect_e8_n105_certificate as prev
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    COS_LOWER_COEFFS,
    GAP_LOWER_COEFFS,
    GAP_UPPER_COEFFS,
    poly_divmod,
    poly_eval,
    positive_no_roots,
)


N = 93
P_NEG = 12
CHEB_DEGREE = 27
MOMENT_MAX = P_NEG + 2 * CHEB_DEGREE + 4
DELTA = Fraction(21, 5)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
S_END = Fraction(248, 5)
T_START = Fraction(30)
T_END = Fraction(969, 20)
EXPECTED_CELLS = 53822
GAP_LOWER_WIDE_COEFFS = (
    GAP_LOWER_COEFFS[0],
    GAP_LOWER_COEFFS[1],
    Fraction(7, 4000),
)
SINE_QUARTIC_WIDE = Fraction(1, 2000)
CHEBYSHEV_INTEGRAL = prev95.CHEBYSHEV_INTEGRAL


def configure() -> None:
    base.N = N
    base.S_START = S_START
    base.S_END = S_END
    base.T_START = T_START
    base.T_END = T_END
    base.GRID_STEP = GRID_STEP


def cap_value() -> Fraction:
    return DELTA * DELTA * 30 * N / (4 * base.D)


def scalar_polynomial_checks() -> dict[str, bool]:
    interval_end = DELTA * DELTA
    gap_lower_residual = [
        Fraction(0),
        Fraction(1) - GAP_LOWER_WIDE_COEFFS[0],
        Fraction(-1, 12) - GAP_LOWER_WIDE_COEFFS[1],
        Fraction(1, 360) - GAP_LOWER_WIDE_COEFFS[2],
        Fraction(-1, 20160),
        Fraction(1, 1814400),
        Fraction(-1, 239500800),
    ]
    gap_upper_residual = [
        Fraction(0),
        GAP_UPPER_COEFFS[0] - 1,
        GAP_UPPER_COEFFS[1] + Fraction(1, 12),
        GAP_UPPER_COEFFS[2] - Fraction(1, 360),
        Fraction(1, 20160),
        Fraction(-1, 1814400),
    ]
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
        "wide_gap_lower_sturm": positive_no_roots(
            gap_lower_residual[1:], Fraction(0), interval_end
        ),
        "wide_gap_upper_sturm": positive_no_roots(
            gap_upper_residual[2:], Fraction(0), interval_end
        ),
        "wide_cos_lower_sturm": positive_no_roots(cos_lower_residual, Fraction(0), interval_end),
    }


def sine_quartic_check() -> bool:
    y_max = DELTA * DELTA / 4
    exponent_poly = [Fraction(0), Fraction(1, 6), Fraction(1, 125)]
    sinc_lower = [Fraction((-1) ** k, factorial(2 * k + 1)) for k in range(8)]
    exp_poly: list[Fraction] = []
    for k in range(6):
        exp_poly = base.poly_add(
            exp_poly,
            base.poly_scale(base.poly_pow(exponent_poly, k), Fraction(1, factorial(k))),
        )
    residual = base.poly_add(base.poly_mul(sinc_lower, exp_poly), [Fraction(-1)])
    quotient = residual[:]
    while quotient and quotient[0] == 0:
        quotient = quotient[1:]
    return positive_no_roots(quotient, Fraction(0), y_max)


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
    )


def moment_integral(power: int) -> int:
    return prev95.moment_integral(power)


def chebyshev_factor() -> list[Fraction]:
    return prev95.chebyshev_factor()


def chebyshev_negative_interval_check() -> bool:
    s_poly = chebyshev_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = prev.helper.poly_mul(x, x)
    xp = prev.helper.poly_pow(x, P_NEG)
    c_neg_x = prev.compose_linear(s_poly, Fraction(-1), Fraction(0))
    lhs = prev.helper.poly_scale(
        prev.helper.poly_mul(
            prev.helper.poly_mul(xp, prev.helper.poly_add(x2, [Fraction(4)])),
            prev.helper.poly_mul(c_neg_x, c_neg_x),
        ),
        Fraction(63, 65) * 16 ** (N - P_NEG),
    )
    rhs = prev.helper.poly_mul(prev.helper.poly_pow(x, N), prev.helper.poly_add(x2, [Fraction(-4)]))
    quotient, remainder = poly_divmod(prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1)), xp)
    quotient2, remainder2 = poly_divmod(quotient, [Fraction(16), Fraction(-1)])
    return (
        not remainder
        and not remainder2
        and poly_eval(quotient2, Fraction(2)) > 0
        and poly_eval(quotient2, Fraction(16)) > 0
        and shared.bernstein_positive_on(quotient2, Fraction(2), Fraction(16))
    )


def chebyshev_positive_interval_check() -> bool:
    s_poly = chebyshev_factor()
    s = [Fraction(0), Fraction(1)]
    s2 = prev.helper.poly_mul(s, s)
    sp = prev.helper.poly_pow(s, P_NEG)
    lhs = prev.helper.poly_scale(
        prev.helper.poly_mul(
            prev.helper.poly_mul(sp, prev.helper.poly_add(s2, [Fraction(4)])),
            prev.helper.poly_mul(s_poly, s_poly),
        ),
        Fraction(63, 65) * 16 ** (N - P_NEG),
    )
    rhs = prev.helper.poly_mul(
        prev.helper.poly_pow(s, N),
        prev.helper.poly_add([Fraction(4)], prev.helper.poly_scale(s2, -1)),
    )
    quotient, remainder = poly_divmod(prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1)), sp)
    return (
        not remainder
        and poly_eval(quotient, Fraction(0)) > 0
        and poly_eval(quotient, Fraction(2)) > 0
        and shared.bernstein_positive_on(quotient, Fraction(0), Fraction(2))
    )


def chebyshev_moment_check() -> bool:
    return (
        moment_integral(12) == 5428655636
        and moment_integral(14) == 1070740575644
        and prev95.chebyshev_integral() == CHEBYSHEV_INTEGRAL
    )


def negative_bound() -> Fraction:
    return Fraction(63, 65) * 16 ** (N - P_NEG) * CHEBYSHEV_INTEGRAL


def old_negative_bound() -> int:
    return (moment_integral(14) + 4 * moment_integral(12)) * 16 ** (N - 12)


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    lower = (
        GAP_LOWER_WIDE_COEFFS[0] * s0
        + GAP_LOWER_WIDE_COEFFS[1] * Fraction(2 * base.D, base.QUARTIC_DENOM) * s0 * s0 / N
        + GAP_LOWER_WIDE_COEFFS[2]
        * Fraction((2 * base.D) ** 2, base.SEXTIC_DENOM)
        * s0**3
        / (N * N)
    )
    upper = (
        GAP_UPPER_COEFFS[0] * t1
        + GAP_UPPER_COEFFS[1] * Fraction(2 * base.D, base.QUARTIC_DENOM) * t1 * t1 / N
        + GAP_UPPER_COEFFS[2]
        * Fraction((2 * base.D) ** 2, base.SEXTIC_DENOM)
        * t1**3
        / (N * N)
    )
    return lower - upper


def local_h_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    return (
        (16 + 240 * COS_LOWER_COEFFS[0]) / (2 * base.D)
        + COS_LOWER_COEFFS[1] * (s1 + t1) / N
        + COS_LOWER_COEFFS[2]
        * Fraction(2 * base.D, base.QUARTIC_DENOM)
        * (s0 * s0 + t0 * t0)
        / (N * N)
        + COS_LOWER_COEFFS[3]
        * Fraction((2 * base.D) ** 2, base.SEXTIC_DENOM)
        * (s1**3 + t1**3)
        / (N**3)
    )


def local_direct_factor_lower(h: Fraction) -> Fraction:
    return Fraction((2 * base.D) ** 2) * h * h - 4


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


def rectangle_fraction_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction
) -> Fraction:
    shape = base.SHAPE
    s_integral = (s1**shape - s0**shape) / shape
    t_integral = (t1**shape - t0**shape) / shape
    return (
        gap_lower(s0, s1, t0, t1) ** 2
        * base.exp_fractional_septic_lower(s1 + t1)
        * s_integral
        * t_integral
        * base.gamma_inv_square()
        / shape
    )


def sine_density_lower(s1: Fraction, t1: Fraction) -> Fraction:
    exponent = (
        2 * base.SINE_QUADRATIC * Fraction(2 * base.D) * (s1 + t1) / N
        + 2
        * SINE_QUARTIC_WIDE
        * Fraction((2 * base.D) ** 2, base.QUARTIC_DENOM)
        * (s1 * s1 + t1 * t1)
        / (N * N)
    )
    return base.exp_fractional_septic_lower(exponent)


def main() -> None:
    configure()
    cells = candidate_cells()
    if len(cells) != EXPECTED_CELLS:
        raise AssertionError(f"expected {EXPECTED_CELLS} retained grid cells, got {len(cells)}")

    scalar_checks = scalar_polynomial_checks()
    gaps: list[Fraction] = []
    h_values: list[Fraction] = []
    direct_values: list[Fraction] = []
    rect_fracs: list[Fraction] = []
    positive_terms: list[Fraction] = []
    sine_exponents: list[Fraction] = []
    a0 = base.a0_lower()
    neg = negative_bound()

    for cell in cells:
        s0, s1, t0, t1 = cell
        gap = gap_lower(s0, s1, t0, t1)
        h = local_h_lower(s0, s1, t0, t1)
        direct = local_direct_factor_lower(h)
        rect_frac = rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = (
            2 * base.SINE_QUADRATIC * Fraction(2 * base.D) * (s1 + t1) / N
            + 2
            * SINE_QUARTIC_WIDE
            * Fraction((2 * base.D) ** 2, base.QUARTIC_DENOM)
            * (s1 * s1 + t1 * t1)
            / (N * N)
        )
        sine = base.exp_fractional_septic_lower(sine_exponent)
        term = (
            a0
            * Fraction((2 * base.D) ** N, N**base.ALPHA)
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

    ratio93 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg
    cheb = chebyshev_factor()

    checks = {
        **scalar_checks,
        "sine_quartic_21_5_sturm": sine_quartic_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_oeis_overlap_0_30": moment_source_check(),
        "shifted_chebyshev_normalization": poly_eval(cheb, Fraction(-16)) == 1,
        "shifted_chebyshev_negative_bernstein": chebyshev_negative_interval_check(),
        "shifted_chebyshev_positive_bernstein": chebyshev_positive_interval_check(),
        "shifted_chebyshev_moment_integral": chebyshev_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h93": min(h_values) > 0,
        "positive_direct93": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio93_gt_1": ratio93 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={DELTA} e_upper={base.E_UPPER}")
    print(
        f"grid_step={GRID_STEP} s_start={S_START} s_end={S_END} "
        f"t_start={T_START} t_end={T_END}"
    )
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("wide_gap_lower_coeffs=99119225/100000000,-7802858/100000000,7/4000")
    print("sine_bound=exp(-z^2/24-z^4/2000)")
    print(
        f"shifted_chebyshev_degree={CHEB_DEGREE} shifted_interval=[7,496] "
        f"moment_max={MOMENT_MAX} cheb_at_0={poly_eval(cheb, Fraction(0))}"
    )
    print(f"cheb_integral={CHEBYSHEV_INTEGRAL}")
    print(f"log10_cheb_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=(63/65)*16^81*cheb_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap93_min={min(gaps)}")
    print(f"h93_min={min(h_values)} direct_factor93_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio93_lower={log10_fraction(ratio93):.2f}")
    print(f"log10_best_cell_ratio93={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=93 rectangular certificate failed")


if __name__ == "__main__":
    main()
