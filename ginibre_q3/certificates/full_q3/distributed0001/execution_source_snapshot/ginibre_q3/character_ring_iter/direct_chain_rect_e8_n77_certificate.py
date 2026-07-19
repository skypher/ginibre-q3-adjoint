#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=77.

This reuses the n=79 rectangular machinery with N=77.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_100, checked against the local m_0..m_70 character-ring log.
This permits a degree-37 square majorant at p=22 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n79_certificate as prev79
import direct_chain_rect_e8_n89_certificate as template
from direct_chain_a0_dominance import log10_fraction


N = 77
P_NEG = 22
NEG_DEGREE = 37
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 57817

NEG_CHEB_COEFFS = (
    Fraction(767948896994567958367, 1000000000000000000000000),
    Fraction(647505183721693979, 250000000000000000000),
    Fraction(141944721854516882869, 100000000000000000000000),
    Fraction(-82298769737381284871, 50000000000000000000000),
    Fraction(-487282091846180087351, 500000000000000000000000),
    Fraction(66121447633777465431, 1000000000000000000000000),
    Fraction(-688958065947111267787, 1000000000000000000000000),
    Fraction(-762088804882658270729, 1000000000000000000000000),
    Fraction(-22104081108837495133, 20000000000000000000000),
    Fraction(-373142301269845078629, 1000000000000000000000000),
    Fraction(145668362854994924077, 100000000000000000000000),
    Fraction(247112122692672018833, 1000000000000000000000000),
    Fraction(-105564968192045461947, 100000000000000000000000),
    Fraction(36919878745841262521, 125000000000000000000000),
    Fraction(1535183699533790551, 3906250000000000000000),
    Fraction(57217939347148309753, 250000000000000000000000),
    Fraction(22491744787497434953, 20000000000000000000000),
    Fraction(62167807984109168449, 125000000000000000000000),
    Fraction(312169157243472886667, 500000000000000000000000),
    Fraction(131295615810682496917, 100000000000000000000000),
    Fraction(-20628268875264016433, 20000000000000000000000),
    Fraction(-79791644592310219087, 50000000000000000000000),
    Fraction(76929114084858071789, 50000000000000000000000),
    Fraction(39560387238531351651, 25000000000000000000000),
    Fraction(-463040974953878536839, 1000000000000000000000000),
    Fraction(-137767064567163737011, 10000000000000000000000000),
    Fraction(11685645466526524959, 62500000000000000000000),
    Fraction(-123670901925416765193, 500000000000000000000000),
    Fraction(862357705779767559779, 1000000000000000000000000),
    Fraction(549095123193334875551, 1000000000000000000000000),
    Fraction(-14794337560290218043, 10000000000000000000000),
    Fraction(-49110253996205664003, 50000000000000000000000),
    Fraction(433419891001757003469, 500000000000000000000000),
    Fraction(55777095022268639129, 250000000000000000000000),
    Fraction(-38343734603905589469, 25000000000000000000000),
    Fraction(-82092794788531803081, 50000000000000000000000),
    Fraction(-86030738014177132197, 125000000000000000000000),
    Fraction(-1099505176235626777, 10000000000000000000000),
)

NEGATIVE_INTEGRAL = Fraction(
    275850426105466974599271222814594639961288531779232471683078843699472794505799470453340725385828976397284323226218809282233398263426608337113956109325735037174032212343506831322288989678916089949771483864475633877431448353313163971042082653208715595538053,
    255150006551666447927467479850110945281873342578074965358424802181970416857399428108290876737920256142293025542128089901049033362686608611252875000499499517287358254803471526768983545297079545658946025000000000000000000000000000000000000000000000000,
)


def install_parameters() -> None:
    template.N = N
    template.P_NEG = P_NEG
    template.NEG_DEGREE = NEG_DEGREE
    template.MOMENT_MAX = MOMENT_MAX
    template.DELTA = DELTA
    template.GRID_STEP = GRID_STEP
    template.S_START = S_START
    template.T_START = T_START
    template.EXPECTED_CELLS = EXPECTED_CELLS
    template.MOMENTS_0_72 = prev79.prev81.MOMENTS_0_100
    template.NEG_CHEB_COEFFS = NEG_CHEB_COEFFS


def configure() -> None:
    install_parameters()
    template.configure()


def cap_value() -> Fraction:
    install_parameters()
    return template.cap_value()


def scalar_polynomial_checks() -> dict[str, bool]:
    install_parameters()
    return {
        name.replace("n89", "n77"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    return prev79.moment_source_check()


def moment_integral(power: int) -> int:
    install_parameters()
    return template.moment_integral(power)


def negative_integral() -> Fraction:
    install_parameters()
    return template.negative_integral()


def negative_bound() -> Fraction:
    install_parameters()
    return template.negative_bound()


def optimized_negative_moment_check() -> bool:
    return (
        moment_integral(22) == 19691620076032366262620
        and moment_integral(24) == 12138032552634345425326212
        and negative_integral() == NEGATIVE_INTEGRAL
    )


def main() -> None:
    configure()
    cells = template.candidate_cells()
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
    a0 = template.base.a0_lower()
    common = Fraction((2 * template.base.D) ** N, N**template.base.ALPHA)
    neg = negative_bound()

    for cell in cells:
        s0, s1, t0, t1 = cell
        gap = template.gap_lower(s0, s1, t0, t1)
        h = template.local_h_lower(s0, s1, t0, t1)
        direct = template.local_direct_factor_lower(h)
        rect_frac = template.rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = template.sine_exponent_upper(s1, t1)
        sine = template.base.exp_fractional_septic_lower(sine_exponent)
        term = a0 * common * rect_frac * h**N * direct * sine
        positive_total += term
        if term > best_term:
            best_term = term
        gaps.append(gap)
        h_values.append(h)
        direct_values.append(direct)
        rect_fracs.append(rect_frac)
        sine_exponents.append(sine_exponent)

    ratio77 = positive_total / neg
    checks = {
        **scalar_checks,
        "sine_sextic_14_3_sturm": template.sine_sextic_check(),
        "root_identity_constants": template.base.root_identity_check(),
        "e_upper_1457_536": template.base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_arxiv_extension": moment_source_check(),
        "optimized_negative_interval_bernstein": template.optimized_negative_interval_check(),
        "optimized_positive_interval_bernstein": template.optimized_positive_interval_check(),
        "optimized_negative_moment_integral": optimized_negative_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h77": min(h_values) > 0,
        "positive_direct77": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio77_gt_1": ratio77 > 1,
    }

    print(f"group={template.base.GROUP} n={N} delta={DELTA} e_upper={template.base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_moment_weight_p={P_NEG} moment_max={MOMENT_MAX}")
    print(f"negative_integral={negative_integral()}")
    print(f"negative_bound=16^{N - P_NEG}*optimized_degree{NEG_DEGREE}_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap77_min={min(gaps)}")
    print(f"h77_min={min(h_values)} direct_factor77_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive77_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative77_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio77_lower={log10_fraction(ratio77):.2f}")
    print(f"log10_best_cell_ratio77={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=77 rectangular certificate failed")


if __name__ == "__main__":
    main()
