#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=105.

This replays the n=113 rectangular machinery with N=105.  The positive
region is a 1/20 by 1/20 rectangular union inside |alpha.v| <= 4.  The
negative-region estimate uses the Chebyshev degree-9 majorant

    (63/65) * 16^(n-8) * S^8 * (S^2+4) * C_9(S)^2,
    C_9(S) = T_9(S/248 - 1) / T_9(-33/31),

for S = X+Y on the negative-sign intervals S in [-16,-2] and S in [0,2].
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n107_certificate as helper
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    poly_divmod,
    poly_eval,
    positive_no_roots,
)


N = 105
P_NEG = 8
CHEB_DEGREE = 9
S_END = Fraction(254, 5)
EXPECTED_CELLS = 63660
CHEBYSHEV_INTEGRAL = Fraction(
    1374322971742349851054687958226748195749248467,
    7483467728620557379409154785380638130176,
)


def configure() -> None:
    base.N = N
    base.S_START = Fraction(40)
    base.S_END = S_END
    base.T_START = Fraction(30)
    base.T_END = Fraction(48)
    base.GRID_STEP = Fraction(1, 20)


def chebyshev_polynomial(degree: int) -> list[Fraction]:
    if degree == 0:
        return [Fraction(1)]
    if degree == 1:
        return [Fraction(0), Fraction(1)]
    previous = [Fraction(1)]
    current = [Fraction(0), Fraction(1)]
    variable = [Fraction(0), Fraction(1)]
    for _ in range(2, degree + 1):
        previous, current = current, helper.poly_add(
            helper.poly_scale(helper.poly_mul(variable, current), 2),
            helper.poly_scale(previous, -1),
        )
    return current


def compose_linear(poly: list[Fraction], scale: Fraction, shift: Fraction) -> list[Fraction]:
    result: list[Fraction] = []
    linear = [shift, scale]
    power = [Fraction(1)]
    for coeff in poly:
        result = helper.poly_add(result, helper.poly_scale(power, coeff))
        power = helper.poly_mul(power, linear)
    return result


def chebyshev_factor() -> list[Fraction]:
    raw = compose_linear(chebyshev_polynomial(CHEB_DEGREE), Fraction(1, 248), Fraction(-1))
    value_at_minus_16 = poly_eval(raw, Fraction(-16))
    return [coeff / value_at_minus_16 for coeff in raw]


def chebyshev_negative_interval_check() -> bool:
    s_poly = chebyshev_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = helper.poly_mul(x, x)
    xp = helper.poly_pow(x, P_NEG)
    c_neg_x = compose_linear(s_poly, Fraction(-1), Fraction(0))
    lhs = helper.poly_scale(
        helper.poly_mul(
            helper.poly_mul(xp, helper.poly_add(x2, [Fraction(4)])),
            helper.poly_mul(c_neg_x, c_neg_x),
        ),
        Fraction(63, 65) * 16 ** (N - P_NEG),
    )
    rhs = helper.poly_mul(helper.poly_pow(x, N), helper.poly_add(x2, [Fraction(-4)]))
    quotient, remainder = poly_divmod(helper.poly_add(lhs, helper.poly_scale(rhs, -1)), xp)
    quotient2, remainder2 = poly_divmod(quotient, [Fraction(16), Fraction(-1)])
    return (
        not remainder
        and not remainder2
        and poly_eval(quotient2, Fraction(2)) > 0
        and poly_eval(quotient2, Fraction(16)) > 0
        and positive_no_roots(quotient2, Fraction(2), Fraction(16))
    )


def chebyshev_positive_interval_check() -> bool:
    s_poly = chebyshev_factor()
    s = [Fraction(0), Fraction(1)]
    s2 = helper.poly_mul(s, s)
    sp = helper.poly_pow(s, P_NEG)
    lhs = helper.poly_scale(
        helper.poly_mul(
            helper.poly_mul(sp, helper.poly_add(s2, [Fraction(4)])),
            helper.poly_mul(s_poly, s_poly),
        ),
        Fraction(63, 65) * 16 ** (N - P_NEG),
    )
    rhs = helper.poly_mul(
        helper.poly_pow(s, N),
        helper.poly_add([Fraction(4)], helper.poly_scale(s2, -1)),
    )
    quotient, remainder = poly_divmod(helper.poly_add(lhs, helper.poly_scale(rhs, -1)), sp)
    return (
        not remainder
        and poly_eval(quotient, Fraction(0)) > 0
        and poly_eval(quotient, Fraction(2)) > 0
        and positive_no_roots(quotient, Fraction(0), Fraction(2))
    )


def chebyshev_integral() -> Fraction:
    factor_square = helper.poly_mul(chebyshev_factor(), chebyshev_factor())
    total = Fraction(0)
    for i, coeff in enumerate(factor_square):
        total += coeff * (
            helper.moment_integral(P_NEG + 2 + i)
            + 4 * helper.moment_integral(P_NEG + i)
        )
    return total


def chebyshev_moment_check() -> bool:
    return (
        helper.moment_integral(12) == 5428655636
        and helper.moment_integral(14) == 1070740575644
        and chebyshev_integral() == CHEBYSHEV_INTEGRAL
    )


def negative_bound() -> Fraction:
    return Fraction(63, 65) * 16 ** (N - P_NEG) * CHEBYSHEV_INTEGRAL


def old_negative_bound() -> int:
    return (helper.moment_integral(14) + 4 * helper.moment_integral(12)) * 16 ** (N - 12)


def main() -> None:
    configure()
    cells = base.candidate_cells()
    if len(cells) != EXPECTED_CELLS:
        raise AssertionError(f"expected {EXPECTED_CELLS} retained grid cells, got {len(cells)}")

    cubic_checks = base.scalar_polynomial_checks()
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
        gap = base.gap_lower(s0, s1, t0, t1)
        h = base.local_h_lower(s0, s1, t0, t1)
        direct = base.local_direct_factor_lower(h)
        rect_frac = base.rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = (
            2 * base.SINE_QUADRATIC * Fraction(2 * base.D) * (s1 + t1) / N
            + 2
            * base.SINE_QUARTIC
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

    ratio105 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg
    cheb = chebyshev_factor()

    checks = {
        **cubic_checks,
        "sine_quartic_sturm": base.sine_scalar_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        "moment_sources_0_30": helper.moment_source_check(),
        "chebyshev_normalization": poly_eval(cheb, Fraction(-16)) == 1,
        "chebyshev_negative_interval_sturm": chebyshev_negative_interval_check(),
        "chebyshev_positive_interval_sturm": chebyshev_positive_interval_check(),
        "chebyshev_moment_integral": chebyshev_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(
            s1 < base.cap_value() and t1 < base.cap_value() for _, s1, _, t1 in cells
        ),
        "positive_gaps": min(gaps) > 0,
        "positive_h105": min(h_values) > 0,
        "positive_direct105": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio105_gt_1": ratio105 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={base.DELTA} e_upper={base.E_UPPER}")
    print(
        f"grid_step={base.GRID_STEP} s_start={base.S_START} s_end={base.S_END} "
        f"t_start={base.T_START} t_end={base.T_END}"
    )
    print(f"retained_cells={len(cells)} cap={base.cap_value()}")
    print(f"chebyshev_degree={CHEB_DEGREE} cheb_at_0={poly_eval(cheb, Fraction(0))}")
    print(f"cheb_integral={CHEBYSHEV_INTEGRAL}")
    print(f"log10_cheb_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=(63/65)*16^97*cheb_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap105_min={min(gaps)}")
    print(f"h105_min={min(h_values)} direct_factor105_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio105_lower={log10_fraction(ratio105):.2f}")
    print(f"log10_best_cell_ratio105={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=105 rectangular certificate failed")


if __name__ == "__main__":
    main()
