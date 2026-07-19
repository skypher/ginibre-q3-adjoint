#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=95.

This replays the n=99,n=97 rectangular machinery with N=95.  The positive
region is a 1/20 by 1/20 rectangular union inside |alpha.v| <= 4.  The
negative-region estimate uses the shifted Chebyshev degree-27 majorant

    (63/65) * 16^(n-12) * S^12 * (S^2+4) * C_27(S)^2,
    C_27(S) = T_27((2S-503)/489) / T_27(-535/489),

for S = X+Y on the negative-sign intervals S in [-16,-2] and S in [0,2].
The interval residuals are checked by exact Bernstein-coefficient positivity
after clearing the non-negative endpoint factors.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import direct_chain_rect_e8_n99_n97_certificate as shared
import direct_chain_rect_e8_n105_certificate as prev
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import poly_divmod, poly_eval


N = 95
P_NEG = 12
CHEB_DEGREE = 27
MOMENT_MAX = P_NEG + 2 * CHEB_DEGREE + 4
S_END = Fraction(919, 20)
T_END = Fraction(907, 20)
EXPECTED_CELLS = 29784
CHEBYSHEV_INTEGRAL = Fraction(
    13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828,
    51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025,
)
MOMENTS_0_70 = shared.MOMENTS_0_64 + (
    180309374169509731329537174820877936398738849331783992400559108044186640253265645187929984,
    12222417127232147684891909975421243170868532113080014550245505153869879217379338376201669320,
    839891813826454630149721108809435050921968630106725839071576059440743027805579950368344519160,
    58491310928223900258448684170711988581617710329632790451360061550785094237630379693882182595544,
    4127047460735629536098477386751185211158919568718740474731125115507267713639183155142366801884608,
    294951001191385645784272903749956304605812964696351517993825555849323854259148399406725458201085948,
)


def configure() -> None:
    base.N = N
    base.S_START = Fraction(40)
    base.S_END = S_END
    base.T_START = Fraction(30)
    base.T_END = T_END
    base.GRID_STEP = Fraction(1, 20)


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(MOMENTS_0_70))
        == MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == MOMENTS_0_70[:31]
    )


def moment_integral(power: int) -> int:
    total = 0
    for k in range(power + 1):
        coeff = comb(power, k)
        total += coeff * MOMENTS_0_70[k + 2] * MOMENTS_0_70[power - k]
        total += coeff * MOMENTS_0_70[k] * MOMENTS_0_70[power - k + 2]
        total -= 2 * coeff * MOMENTS_0_70[k + 1] * MOMENTS_0_70[power - k + 1]
    return total


def chebyshev_factor() -> list[Fraction]:
    raw = prev.compose_linear(
        prev.chebyshev_polynomial(CHEB_DEGREE),
        Fraction(2, 489),
        Fraction(-503, 489),
    )
    value_at_minus_16 = poly_eval(raw, Fraction(-16))
    return [coeff / value_at_minus_16 for coeff in raw]


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


def chebyshev_integral() -> Fraction:
    factor_square = prev.helper.poly_mul(chebyshev_factor(), chebyshev_factor())
    total = Fraction(0)
    for i, coeff in enumerate(factor_square):
        total += coeff * (
            moment_integral(P_NEG + 2 + i)
            + 4 * moment_integral(P_NEG + i)
        )
    return total


def chebyshev_moment_check() -> bool:
    return (
        moment_integral(12) == 5428655636
        and moment_integral(14) == 1070740575644
        and chebyshev_integral() == CHEBYSHEV_INTEGRAL
    )


def negative_bound() -> Fraction:
    return Fraction(63, 65) * 16 ** (N - P_NEG) * CHEBYSHEV_INTEGRAL


def old_negative_bound() -> int:
    return (moment_integral(14) + 4 * moment_integral(12)) * 16 ** (N - 12)


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

    ratio95 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg
    cheb = chebyshev_factor()

    checks = {
        **cubic_checks,
        "sine_quartic_sturm": base.sine_scalar_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_oeis_overlap_0_30": moment_source_check(),
        "shifted_chebyshev_normalization": poly_eval(cheb, Fraction(-16)) == 1,
        "shifted_chebyshev_negative_bernstein": chebyshev_negative_interval_check(),
        "shifted_chebyshev_positive_bernstein": chebyshev_positive_interval_check(),
        "shifted_chebyshev_moment_integral": chebyshev_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(
            s1 < base.cap_value() and t1 < base.cap_value() for _, s1, _, t1 in cells
        ),
        "positive_gaps": min(gaps) > 0,
        "positive_h95": min(h_values) > 0,
        "positive_direct95": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio95_gt_1": ratio95 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={base.DELTA} e_upper={base.E_UPPER}")
    print(
        f"grid_step={base.GRID_STEP} s_start={base.S_START} s_end={base.S_END} "
        f"t_start={base.T_START} t_end={base.T_END}"
    )
    print(f"retained_cells={len(cells)} cap={base.cap_value()}")
    print(
        f"shifted_chebyshev_degree={CHEB_DEGREE} shifted_interval=[7,496] "
        f"moment_max={MOMENT_MAX} cheb_at_0={poly_eval(cheb, Fraction(0))}"
    )
    print(f"cheb_integral={CHEBYSHEV_INTEGRAL}")
    print(f"log10_cheb_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=(63/65)*16^83*cheb_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap95_min={min(gaps)}")
    print(f"h95_min={min(h_values)} direct_factor95_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio95_lower={log10_fraction(ratio95):.2f}")
    print(f"log10_best_cell_ratio95={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=95 rectangular certificate failed")


if __name__ == "__main__":
    main()
