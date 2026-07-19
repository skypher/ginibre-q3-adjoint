#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=89.

The positive side uses the n=91 rectangular machinery with root-angle cap
|alpha.v| <= 14/3, grid step 1/20, a cap-adjusted cubic gap coefficient,
and the Weyl-sine scalar bound

    log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2880 - z^6/103800.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,m_72, checked against the local m_0..m_70 character-ring log.
This permits a degree-28 square majorant at the same p=12 moment weight.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial

import direct_chain_rect_e8_n99_n97_certificate as shared
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n105_certificate as prev
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    COS_LOWER_COEFFS,
    GAP_UPPER_COEFFS,
    poly_eval,
    positive_no_roots,
)


N = 89
P_NEG = 12
NEG_DEGREE = 28
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 130740
GAP_LOWER_N89_COEFFS = (
    Fraction(99119225, 100000000),
    Fraction(-7802858, 100000000),
    Fraction(169078, 100000000),
)
SINE_QUARTIC_N89 = Fraction(1, 2880)
SINE_SEXTIC_N89 = Fraction(1, 103800)
MOMENTS_0_72 = prev95.MOMENTS_0_70 + (
    21345605852669727193269624795683377626642321884441139772528810920344326036087847049938869351873014236,
    1563884705668584068925407561523930788789172461189673466806410484540072278551845672660906888535095697136,
)

NEG_CHEB_COEFFS = (
    Fraction("-0.007429629917251921"),
    Fraction("0.005788235384571325"),
    Fraction("-0.0007816247669013144"),
    Fraction("-0.007558047667202726"),
    Fraction("0.017058885120560974"),
    Fraction("-0.013658151971889694"),
    Fraction("-0.014841460395873158"),
    Fraction("0.020325497785555396"),
    Fraction("-0.007127888530601548"),
    Fraction("-0.00930930288786367"),
    Fraction("0.012553954606718088"),
    Fraction("-0.006511367168472324"),
    Fraction("-0.0025894374850525116"),
    Fraction("0.006882175720294705"),
    Fraction("-0.00491803277544681"),
    Fraction("0.0019916480971363163"),
    Fraction("0.0036180996067065666"),
    Fraction("-0.00338935244934804"),
    Fraction("0.0020052097887612707"),
    Fraction("0.0005084796355251987"),
    Fraction("-0.0022611393291781354"),
    Fraction("0.0016868903540933997"),
    Fraction("0.00023754349461354313"),
    Fraction("-0.0008099507722295414"),
    Fraction("0.0009292290349050409"),
    Fraction("-8.235052353854276e-05"),
    Fraction("-0.0004559837369155558"),
    Fraction("0.00043121022673966845"),
    Fraction("0.0002959848125583796"),
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
        Fraction(1) - GAP_LOWER_N89_COEFFS[0],
        Fraction(-1, 12) - GAP_LOWER_N89_COEFFS[1],
        Fraction(1, 360) - GAP_LOWER_N89_COEFFS[2],
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
        "n89_gap_lower_sturm": positive_no_roots(
            gap_lower_residual[1:], Fraction(0), interval_end
        ),
        "n89_gap_upper_sturm": positive_no_roots(
            gap_upper_residual[2:], Fraction(0), interval_end
        ),
        "n89_cos_lower_sturm": positive_no_roots(
            cos_lower_residual, Fraction(0), interval_end
        ),
    }


def sine_sextic_check() -> bool:
    y_max = DELTA * DELTA / 4
    exponent_poly = [Fraction(0), Fraction(1, 6), Fraction(1, 180), 64 * SINE_SEXTIC_N89]
    sinc_lower = [Fraction((-1) ** k, factorial(2 * k + 1)) for k in range(8)]
    exp_poly: list[Fraction] = []
    for k in range(8):
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
    arxiv_moments: dict[int, int] = {}
    with (here.parent / "references" / "arxiv_2412_21189_e8_m71_m72.txt").open() as handle:
        for line in handle:
            match = shared.re.search(r"\bm_\s*(\d+)\s*=\s+(-?\d+)", line)
            if match:
                arxiv_moments[int(match.group(1))] = int(match.group(2))
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
        and arxiv_moments == {71: MOMENTS_0_72[71], 72: MOMENTS_0_72[72]}
    )


def moment_integral(power: int) -> int:
    total = 0
    for k in range(power + 1):
        coeff = comb(power, k)
        total += coeff * MOMENTS_0_72[k + 2] * MOMENTS_0_72[power - k]
        total += coeff * MOMENTS_0_72[k] * MOMENTS_0_72[power - k + 2]
        total -= 2 * coeff * MOMENTS_0_72[k + 1] * MOMENTS_0_72[power - k + 1]
    return total


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
            110574439399728242778064347972737777648628809018978860670416750279550250166056988167096970632724370157057840918236970893796884473467113003625506475576422917439930098480783344439652958219223,
            62390165395220237149681115433261077065698483356952027438641538929482648026739857327216912260160826558355727591021787674134419528043283644560507030415015625000000000000000000000000000000000,
        )
    )


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    lower = (
        GAP_LOWER_N89_COEFFS[0] * s0
        + GAP_LOWER_N89_COEFFS[1] * Fraction(2 * base.D, base.QUARTIC_DENOM) * s0 * s0 / N
        + GAP_LOWER_N89_COEFFS[2]
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
        * SINE_QUARTIC_N89
        * Fraction((2 * base.D) ** 2, base.QUARTIC_DENOM)
        * (s1 * s1 + t1 * t1)
        / (N * N)
        + 2
        * SINE_SEXTIC_N89
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

    ratio89 = positive_total / neg
    checks = {
        **scalar_checks,
        "sine_sextic_14_3_sturm": sine_sextic_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_arxiv_extension": moment_source_check(),
        "optimized_negative_interval_bernstein": optimized_negative_interval_check(),
        "optimized_positive_interval_bernstein": optimized_positive_interval_check(),
        "optimized_negative_moment_integral": optimized_negative_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h89": min(h_values) > 0,
        "positive_direct89": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio89_gt_1": ratio89 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={DELTA} e_upper={base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^77*optimized_degree28_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap89_min={min(gaps)}")
    print(f"h89_min={min(h_values)} direct_factor89_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive89_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative89_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio89_lower={log10_fraction(ratio89):.2f}")
    print(f"log10_best_cell_ratio89={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=89 rectangular certificate failed")


if __name__ == "__main__":
    main()
