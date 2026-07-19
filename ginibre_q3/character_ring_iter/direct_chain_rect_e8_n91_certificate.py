#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=91.

The positive side reuses the n=93 rectangular machinery with root-angle cap
|alpha.v| <= 13/3 and a sharper Weyl-sine scalar bound

    log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2880 - z^6/115000.

The negative side keeps the degree-27, moment-70 ceiling, but replaces the
single Chebyshev minimax factor by a weighted moment-norm square majorant.
The two majorant inequalities are checked by exact Bernstein coefficients.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import direct_chain_rect_e8_n99_n97_certificate as shared
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n105_certificate as prev
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    COS_LOWER_COEFFS,
    GAP_LOWER_COEFFS,
    GAP_UPPER_COEFFS,
    poly_eval,
    positive_no_roots,
)


N = 91
P_NEG = 12
NEG_DEGREE = 27
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(13, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 69502
GAP_LOWER_N91_COEFFS = (
    GAP_LOWER_COEFFS[0],
    GAP_LOWER_COEFFS[1],
    Fraction(1562, 890661),
)
SINE_QUARTIC_N91 = Fraction(1, 2880)
SINE_SEXTIC_N91 = Fraction(1, 115000)

NEG_CHEB_COEFFS = (
    Fraction("-5.7967396718931195e-06") + Fraction(1, 100000000),
    Fraction("5.174666956979492e-06"),
    Fraction("-3.4425906740311743e-06"),
    Fraction("9.758209880849629e-07"),
    Fraction("1.685371993460304e-06"),
    Fraction("-3.94656659700228e-06"),
    Fraction("5.283353156116014e-06"),
    Fraction("-5.354690976028374e-06"),
    Fraction("4.082781085947319e-06"),
    Fraction("-1.6824617082595656e-06"),
    Fraction("-1.367920359192741e-06"),
    Fraction("4.411521214126105e-06"),
    Fraction("-6.740666666943214e-06"),
    Fraction("7.745312968246953e-06"),
    Fraction("-7.055293622565833e-06"),
    Fraction("4.6486310655844e-06"),
    Fraction("-9.020868778249604e-07"),
    Fraction("-3.432024191780943e-06"),
    Fraction("7.328326454174731e-06"),
    Fraction("-9.664759274955677e-06"),
    Fraction("9.45809259108623e-06"),
    Fraction("-6.138448573892445e-06"),
    Fraction("-1.567429287892005e-07"),
    Fraction("8.280781035226091e-06"),
    Fraction("-1.5762129488499434e-05"),
    Fraction("1.8442569503947575e-05"),
    Fraction("-9.997590178203355e-06"),
    Fraction("-3.64528811441427e-05"),
)


def configure() -> None:
    base.N = N
    base.S_START = S_START
    base.S_END = cap_value()
    base.T_START = T_START
    base.T_END = cap_value()
    base.GRID_STEP = GRID_STEP


def cap_value() -> Fraction:
    return DELTA * DELTA * 30 * N / (4 * base.D)


def scalar_polynomial_checks() -> dict[str, bool]:
    interval_end = DELTA * DELTA
    gap_lower_residual = [
        Fraction(0),
        Fraction(1) - GAP_LOWER_N91_COEFFS[0],
        Fraction(-1, 12) - GAP_LOWER_N91_COEFFS[1],
        Fraction(1, 360) - GAP_LOWER_N91_COEFFS[2],
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
        "n91_gap_lower_sturm": positive_no_roots(
            gap_lower_residual[1:], Fraction(0), interval_end
        ),
        "n91_gap_upper_sturm": positive_no_roots(
            gap_upper_residual[2:], Fraction(0), interval_end
        ),
        "n91_cos_lower_sturm": positive_no_roots(
            cos_lower_residual, Fraction(0), interval_end
        ),
    }


def sine_sextic_check() -> bool:
    y_max = DELTA * DELTA / 4
    exponent_poly = [Fraction(0), Fraction(1, 6), Fraction(1, 180), Fraction(16, 28750)]
    sinc_lower = [Fraction((-1) ** k, factorial(2 * k + 1)) for k in range(8)]
    exp_poly: list[Fraction] = []
    for k in range(7):
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


def negative_factor() -> list[Fraction]:
    poly: list[Fraction] = []
    for degree, coeff in enumerate(NEG_CHEB_COEFFS):
        basis = prev.compose_linear(
            prev.chebyshev_polynomial(degree),
            Fraction(2, 489),
            Fraction(-503, 489),
        )
        poly = prev.helper.poly_add(poly, prev.helper.poly_scale(basis, coeff))
    return poly


def negative_integral() -> Fraction:
    factor_square = prev.helper.poly_mul(negative_factor(), negative_factor())
    total = Fraction(0)
    for i, coeff in enumerate(factor_square):
        total += coeff * (
            moment_integral(P_NEG + 2 + i) + 4 * moment_integral(P_NEG + i)
        )
    return total


def negative_bound() -> Fraction:
    return 16 ** (N - P_NEG) * negative_integral()


def old_negative_bound() -> int:
    return (moment_integral(14) + 4 * moment_integral(12)) * 16 ** (N - 12)


def optimized_negative_interval_check() -> bool:
    factor = negative_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = prev.helper.poly_mul(x, x)
    p_neg_x = prev.compose_linear(factor, Fraction(-1), Fraction(0))
    lhs = prev.helper.poly_scale(
        prev.helper.poly_mul(
            prev.helper.poly_add(x2, [Fraction(4)]),
            prev.helper.poly_mul(p_neg_x, p_neg_x),
        ),
        16 ** (N - P_NEG),
    )
    rhs = prev.helper.poly_mul(
        prev.helper.poly_pow(x, N - P_NEG),
        prev.helper.poly_add(x2, [Fraction(-4)]),
    )
    residual = prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1))
    return (
        poly_eval(residual, Fraction(2)) > 0
        and poly_eval(residual, Fraction(16)) > 0
        and shared.bernstein_positive_on(residual, Fraction(2), Fraction(16))
    )


def optimized_positive_interval_check() -> bool:
    factor = negative_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = prev.helper.poly_mul(x, x)
    lhs = prev.helper.poly_scale(
        prev.helper.poly_mul(
            prev.helper.poly_add(x2, [Fraction(4)]),
            prev.helper.poly_mul(factor, factor),
        ),
        16 ** (N - P_NEG),
    )
    rhs = prev.helper.poly_mul(
        prev.helper.poly_pow(x, N - P_NEG),
        prev.helper.poly_add([Fraction(4)], prev.helper.poly_scale(x2, -1)),
    )
    residual = prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1))
    return (
        poly_eval(residual, Fraction(0)) > 0
        and poly_eval(residual, Fraction(2)) > 0
        and shared.bernstein_positive_on(residual, Fraction(0), Fraction(2))
    )


def optimized_negative_moment_check() -> bool:
    return (
        moment_integral(12) == 5428655636
        and moment_integral(14) == 1070740575644
        and negative_integral()
        == Fraction(
            154899922705202397910663859891726275260116389741415466199911998543539776953588818472524235571552289212256879365191595587107363092116739932995012713161412017164982537410758765023648674235489,
            3865401811617780810470268062535894321947867605658317875555118084768343566134415279121319626970868866380880938847280490359795809438691655382933018750000000000000000000000000000000000000000,
        )
    )


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    lower = (
        GAP_LOWER_N91_COEFFS[0] * s0
        + GAP_LOWER_N91_COEFFS[1] * Fraction(2 * base.D, base.QUARTIC_DENOM) * s0 * s0 / N
        + GAP_LOWER_N91_COEFFS[2]
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
    while s0 + GRID_STEP <= cap:
        t0 = T_START
        while t0 + GRID_STEP <= cap:
            s1 = s0 + GRID_STEP
            t1 = t0 + GRID_STEP
            gap = gap_lower(s0, s1, t0, t1)
            if gap > 0:
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


def sine_exponent_upper(s1: Fraction, t1: Fraction) -> Fraction:
    return (
        2 * base.SINE_QUADRATIC * Fraction(2 * base.D) * (s1 + t1) / N
        + 2
        * SINE_QUARTIC_N91
        * Fraction((2 * base.D) ** 2, base.QUARTIC_DENOM)
        * (s1 * s1 + t1 * t1)
        / (N * N)
        + 2
        * SINE_SEXTIC_N91
        * Fraction((2 * base.D) ** 2, base.SEXTIC_DENOM)
        * (s1**3 + t1**3)
        / (N**3)
    )


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
    sine_exponents: list[Fraction] = []
    positive_total = Fraction(0)
    best_term = Fraction(0)
    a0 = base.a0_lower()
    common = Fraction((2 * base.D) ** N, N**base.ALPHA)
    neg = negative_bound()

    for cell in cells:
        s0, s1, t0, t1 = cell
        gap = gap_lower(s0, s1, t0, t1)
        h = local_h_lower(s0, s1, t0, t1)
        direct = local_direct_factor_lower(h)
        rect_frac = rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = sine_exponent_upper(s1, t1)
        sine = base.exp_fractional_septic_lower(sine_exponent)
        term = a0 * common * rect_frac * h**N * direct * sine
        positive_total += term
        if term > best_term:
            best_term = term
        gaps.append(gap)
        h_values.append(h)
        direct_values.append(direct)
        rect_fracs.append(rect_frac)
        sine_exponents.append(sine_exponent)

    ratio91 = positive_total / neg
    checks = {
        **scalar_checks,
        "sine_sextic_13_3_sturm": sine_sextic_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_oeis_overlap_0_30": moment_source_check(),
        "optimized_negative_interval_bernstein": optimized_negative_interval_check(),
        "optimized_positive_interval_bernstein": optimized_positive_interval_check(),
        "optimized_negative_moment_integral": optimized_negative_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h91": min(h_values) > 0,
        "positive_direct91": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio91_gt_1": ratio91 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={DELTA} e_upper={base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,1562/890661")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/115000)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^79*optimized_degree27_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap91_min={min(gaps)}")
    print(f"h91_min={min(h_values)} direct_factor91_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive91_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative91_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio91_lower={log10_fraction(ratio91):.2f}")
    print(f"log10_best_cell_ratio91={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=91 rectangular certificate failed")


if __name__ == "__main__":
    main()
