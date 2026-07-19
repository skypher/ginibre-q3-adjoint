#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=79.

This reuses the n=81 rectangular machinery with N=79.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_100, checked against the local m_0..m_70 character-ring log.
This permits a degree-37 square majorant at p=22 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n81_certificate as prev81
import direct_chain_rect_e8_n89_certificate as template
from direct_chain_a0_dominance import log10_fraction


N = 79
P_NEG = 22
NEG_DEGREE = 37
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 68095

NEG_CHEB_COEFFS = (
    Fraction(-31663311737934886999, 50000000000000000000000),
    Fraction(-48577487645705727093, 50000000000000000000000),
    Fraction(-71927003886043674401, 100000000000000000000000),
    Fraction(-22171000412973586867, 50000000000000000000000),
    Fraction(15205976837429048473, 100000000000000000000000),
    Fraction(35951807895168239327, 500000000000000000000000),
    Fraction(-31596187065533944557, 100000000000000000000000),
    Fraction(-19462550499421748059, 100000000000000000000000),
    Fraction(-16937043295424458877, 1000000000000000000000000),
    Fraction(49760448475591411943, 100000000000000000000000),
    Fraction(12574671886098762451, 10000000000000000000000),
    Fraction(20499218709752955667, 20000000000000000000000),
    Fraction(33012766362614695041, 500000000000000000000000),
    Fraction(-399931493631060567, 2000000000000000000000),
    Fraction(36913769116916839631, 100000000000000000000000),
    Fraction(32282739591410566533, 100000000000000000000000),
    Fraction(-10711522480626099913, 20000000000000000000000),
    Fraction(-15798519512098445969, 100000000000000000000000),
    Fraction(87692277458233276741, 100000000000000000000000),
    Fraction(17643859035040676647, 100000000000000000000000),
    Fraction(-33515148572339565527, 50000000000000000000000),
    Fraction(-46665961791752983647, 100000000000000000000000),
    Fraction(-1142486825771279351, 1000000000000000000000),
    Fraction(-14726513193144257857, 10000000000000000000000),
    Fraction(-29963380890479561987, 500000000000000000000000),
    Fraction(3046729976622868921, 10000000000000000000000),
    Fraction(-313403249324376233, 625000000000000000000),
    Fraction(32044149345630445397, 500000000000000000000000),
    Fraction(89464285851346516451, 100000000000000000000000),
    Fraction(57135208592176175549, 100000000000000000000000),
    Fraction(17120821880204376763, 100000000000000000000000),
    Fraction(-26797313343402760609, 100000000000000000000000),
    Fraction(-25306136695862691867, 25000000000000000000000),
    Fraction(-1858351971267914243, 2000000000000000000000),
    Fraction(-5633845538578637149, 2000000000000000000000000),
    Fraction(4790257118526445137, 10000000000000000000000),
    Fraction(28722732222587237439, 100000000000000000000000),
    Fraction(55850628756710893501, 1000000000000000000000000),
)

NEGATIVE_INTEGRAL = Fraction(
    29822698848482859512999327002571175457028680982472229890974727314690982554293289858766526952918578162868771228943226337357855616889719317471018088299258706043797188858166167899916127307423211905059325604991005365894898566191400370099111060755564260674533,
    10206000262066657917098699194004437811274933703122998614336992087278816674295977124331635069516810245691721021685123596041961334507464344450115000019979980691494330192138861070759341811883181826357841000000000000000000000000000000000000000000000000,
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
    template.MOMENTS_0_72 = prev81.MOMENTS_0_100
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
        name.replace("n89", "n79"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    return prev81.moment_source_check()


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

    ratio79 = positive_total / neg
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
        "positive_h79": min(h_values) > 0,
        "positive_direct79": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio79_gt_1": ratio79 > 1,
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
    print(f"gap79_min={min(gaps)}")
    print(f"h79_min={min(h_values)} direct_factor79_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive79_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative79_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio79_lower={log10_fraction(ratio79):.2f}")
    print(f"log10_best_cell_ratio79={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=79 rectangular certificate failed")


if __name__ == "__main__":
    main()
