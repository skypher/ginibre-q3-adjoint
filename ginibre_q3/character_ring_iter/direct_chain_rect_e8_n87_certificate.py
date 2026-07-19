#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=87.

This reuses the n=89 rectangular machinery with N=87.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_76, checked against the local m_0..m_70 character-ring log.
This permits a degree-30 square majorant at the same p=12 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n89_certificate as template
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n99_n97_certificate as shared
from direct_chain_a0_dominance import log10_fraction


N = 87
P_NEG = 12
NEG_DEGREE = 30
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 116526
MOMENTS_0_76 = prev95.MOMENTS_0_70 + (
    21345605852669727193269624795683377626642321884441139772528810920344326036087847049938869351873014236,
    1563884705668584068925407561523930788789172461189673466806410484540072278551845672660906888535095697136,
    115966203895188154500225464991952543976209799949711939008096002629888564857244232623215412510815691511056,
    8701310967044097147019266577057644192169884031190810920786174911051333498536424762493837189576911591909200,
    660486089999020109738643643713150053248192927070524315500997576912261631868550532947681596754927281427810448,
    50707311959735108441198566662019983734164705456998999287374783721659669815172185037522981492125524895220454760,
)

NEG_CHEB_COEFFS = (
    Fraction("0.0016725178130911106"),
    Fraction("-0.00144326947757074386"),
    Fraction("0.000702410045842463714"),
    Fraction("0.000740283426469716638"),
    Fraction("-0.00308009002164543475"),
    Fraction("0.0049893834687314122"),
    Fraction("-0.000692786698166533623"),
    Fraction("-0.00374394886987779998"),
    Fraction("0.00376417124010575447"),
    Fraction("-0.00219557302036715212"),
    Fraction("0.000183996616236802897"),
    Fraction("0.00181274478738766508"),
    Fraction("-0.0031953377141447966"),
    Fraction("0.00305043208672802499"),
    Fraction("-0.00125123414534050975"),
    Fraction("-0.00205164844955836676"),
    Fraction("0.00482910069535850538"),
    Fraction("-0.0038586499472885896"),
    Fraction("-0.000199393106972616793"),
    Fraction("0.00438521803108860107"),
    Fraction("-0.00553747542004831301"),
    Fraction("0.00235639173453197235"),
    Fraction("0.0037272836409665052"),
    Fraction("-0.00755291809341343112"),
    Fraction("0.00337838011008285153"),
    Fraction("0.00611166659834407652"),
    Fraction("-0.00929938615252917135"),
    Fraction("-0.00316945021349966183"),
    Fraction("0.0123798465097602201"),
    Fraction("0.0103377711445294344"),
    Fraction("0.00238718045982216258"),
)

NEGATIVE_INTEGRAL = Fraction(
    93585328496270081696599778390722481293767909668608279570199780293814811208211678612684601640822568478824012644817918965028087964159078576255670743626137222438605197763195972232741016641340705851985699,
    704671224371538857426704268543001678147262865304563156507708291095378414199985724395911100200841907504391374851431646822661155836061171868067289422290307850980250000000000000000000000000000000000000000,
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
    template.MOMENTS_0_72 = MOMENTS_0_76
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
        name.replace("n89", "n87"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    arxiv_moments: dict[int, int] = {}
    with (here.parent / "references" / "arxiv_2412_21189_e8_m71_m76.txt").open() as handle:
        for line in handle:
            match = shared.re.search(r"\bm_\s*(\d+)\s*=\s+(-?\d+)", line)
            if match:
                arxiv_moments[int(match.group(1))] = int(match.group(2))
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
        and arxiv_moments == {k: MOMENTS_0_76[k] for k in range(71, 77)}
    )


def moment_integral(power: int) -> int:
    install_parameters()
    return template.moment_integral(power)


def negative_integral() -> Fraction:
    install_parameters()
    return template.negative_integral()


def negative_bound() -> Fraction:
    install_parameters()
    return template.negative_bound()


def old_negative_bound() -> int:
    install_parameters()
    return template.old_negative_bound()


def optimized_negative_moment_check() -> bool:
    return (
        moment_integral(12) == 5428655636
        and moment_integral(14) == 1070740575644
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

    ratio87 = positive_total / neg
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
        "positive_h87": min(h_values) > 0,
        "positive_direct87": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio87_gt_1": ratio87 > 1,
    }

    print(f"group={template.base.GROUP} n={N} delta={DELTA} e_upper={template.base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^75*optimized_degree30_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap87_min={min(gaps)}")
    print(f"h87_min={min(h_values)} direct_factor87_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive87_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative87_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio87_lower={log10_fraction(ratio87):.2f}")
    print(f"log10_best_cell_ratio87={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=87 rectangular certificate failed")


if __name__ == "__main__":
    main()
