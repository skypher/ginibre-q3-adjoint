#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=85.

This reuses the n=87 rectangular machinery with N=85.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_84, checked against the local m_0..m_70 character-ring log.
This permits a degree-34 square majorant at the same p=12 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n89_certificate as template
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n99_n97_certificate as shared
from direct_chain_a0_dominance import log10_fraction


N = 85
P_NEG = 12
NEG_DEGREE = 34
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 103495
POSITIVE_BERNSTEIN_PIECES = 32
MOMENTS_0_84 = prev95.MOMENTS_0_70 + (
    21345605852669727193269624795683377626642321884441139772528810920344326036087847049938869351873014236,
    1563884705668584068925407561523930788789172461189673466806410484540072278551845672660906888535095697136,
    115966203895188154500225464991952543976209799949711939008096002629888564857244232623215412510815691511056,
    8701310967044097147019266577057644192169884031190810920786174911051333498536424762493837189576911591909200,
    660486089999020109738643643713150053248192927070524315500997576912261631868550532947681596754927281427810448,
    50707311959735108441198566662019983734164705456998999287374783721659669815172185037522981492125524895220454760,
    3936502147198070963491467364963401744755672920873807463653657729434923374395918948001914452306516369775664127648,
    308952142346690731007140204134530998740732421980595804682008261520522367193737777061840984047494093828899384184168,
    24508866165177630891480898882373256902087356919011778926541309325536030302661368517116941510852497647544401414601460,
    1964805416254688176296872219010410027260184853056137299893399904773235151901215515742168232755697937844307110882046648,
    159145978587886029421091500725137512347851291352847196386830845325516385844881545558708812608145902846557201716848974752,
    13021752184433160337313743494828710709292098762448706109681625041927502515163426233897278123992631639608974205538547354348,
    1076120557185963103254375093999197941409443803061870402487867373186878754713606191533570405896237096777471987815243734312122,
    89803336303524881408230011577463211534753145023307206844345162601368469336098768993290975115793996843772571568866707086318102,
)

NEG_CHEB_COEFFS = (
    Fraction("-752720140400008801/1000000000000000000000"),
    Fraction("646318179588940833/1000000000000000000000"),
    Fraction("-303553462921636217/1000000000000000000000"),
    Fraction("-349825856605592183/1000000000000000000000"),
    Fraction("138771476892079787/100000000000000000000"),
    Fraction("-8717769410135693/4000000000000000000"),
    Fraction("48204132982028003/200000000000000000000"),
    Fraction("82221412899908489/50000000000000000000"),
    Fraction("-6388704404615541/4000000000000000000"),
    Fraction("900514315058171437/1000000000000000000000"),
    Fraction("-71557015341184853/1250000000000000000000"),
    Fraction("-45431042962154829/62500000000000000000"),
    Fraction("24917962369789999/20000000000000000000"),
    Fraction("-124674317078151309/100000000000000000000"),
    Fraction("109710343527977261/200000000000000000000"),
    Fraction("665701991850927441/1000000000000000000000"),
    Fraction("-165945389563715209/100000000000000000000"),
    Fraction("79587984221067033/50000000000000000000"),
    Fraction("-46087433666614393/125000000000000000000"),
    Fraction("-59507661115958371/50000000000000000000"),
    Fraction("49957276375438037/25000000000000000000"),
    Fraction("-8893270877328611/6250000000000000000"),
    Fraction("-140585469200787581/500000000000000000000"),
    Fraction("48869650415203159/25000000000000000000"),
    Fraction("-108877673654480503/50000000000000000000"),
    Fraction("437202782416299209/1000000000000000000000"),
    Fraction("196245116411087077/100000000000000000000"),
    Fraction("-32030890287335919/12500000000000000000"),
    Fraction("28389688917118119/100000000000000000000"),
    Fraction("252874506859244183/100000000000000000000"),
    Fraction("-25305375637425017/12500000000000000000"),
    Fraction("-10160378686874213/6250000000000000000"),
    Fraction("123030686054563911/50000000000000000000"),
    Fraction("60588271568714927/25000000000000000000"),
    Fraction("4550586000594209/7812500000000000000"),
)

NEGATIVE_INTEGRAL = Fraction(
    65148687217721097618389819413073654708609879754137122790264611448139731858033630148113111544815479392877378205889026723166072799733022265663747429139658760314718523759506158358410320170203871509898553101174638646758322907219,
    4665330740770604793149491012159339221112632191672044571621781307705608579427809365087496794344569859695546279931419517931942002264051301274362377419277716193612079960316166376689253006250000000000000000000000000000000000000000,
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
    template.MOMENTS_0_72 = MOMENTS_0_84
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
        name.replace("n89", "n85"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    arxiv_moments: dict[int, int] = {}
    with (here.parent / "references" / "arxiv_2412_21189_e8_m71_m84.txt").open() as handle:
        for line in handle:
            match = shared.re.search(r"\bm_\s*(\d+)\s*=\s+(-?\d+)", line)
            if match:
                arxiv_moments[int(match.group(1))] = int(match.group(2))
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
        and arxiv_moments == {k: MOMENTS_0_84[k] for k in range(71, 85)}
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


def optimized_positive_subdivided_interval_check() -> bool:
    install_parameters()
    factor = template.negative_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = template.prev.helper.poly_mul(x, x)
    lhs = template.prev.helper.poly_scale(
        template.prev.helper.poly_mul(
            template.prev.helper.poly_add(x2, [Fraction(4)]),
            template.prev.helper.poly_mul(factor, factor),
        ),
        16 ** (N - P_NEG),
    )
    rhs = template.prev.helper.poly_mul(
        template.prev.helper.poly_pow(x, N - P_NEG),
        template.prev.helper.poly_add([Fraction(4)], template.prev.helper.poly_scale(x2, -1)),
    )
    residual = template.prev.helper.poly_add(lhs, template.prev.helper.poly_scale(rhs, -1))
    width = Fraction(2, POSITIVE_BERNSTEIN_PIECES)
    return (
        template.poly_eval(residual, Fraction(0)) > 0
        and template.poly_eval(residual, Fraction(2)) > 0
        and all(
            shared.bernstein_positive_on(residual, i * width, (i + 1) * width)
            for i in range(POSITIVE_BERNSTEIN_PIECES)
        )
    )


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

    ratio85 = positive_total / neg
    checks = {
        **scalar_checks,
        "sine_sextic_14_3_sturm": template.sine_sextic_check(),
        "root_identity_constants": template.base.root_identity_check(),
        "e_upper_1457_536": template.base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_arxiv_extension": moment_source_check(),
        "optimized_negative_interval_bernstein": template.optimized_negative_interval_check(),
        f"optimized_positive_interval_bernstein_{POSITIVE_BERNSTEIN_PIECES}": (
            optimized_positive_subdivided_interval_check()
        ),
        "optimized_negative_moment_integral": optimized_negative_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h85": min(h_values) > 0,
        "positive_direct85": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio85_gt_1": ratio85 > 1,
    }

    print(f"group={template.base.GROUP} n={N} delta={DELTA} e_upper={template.base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"positive_interval_subdivision_count={POSITIVE_BERNSTEIN_PIECES}")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^73*optimized_degree34_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap85_min={min(gaps)}")
    print(f"h85_min={min(h_values)} direct_factor85_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive85_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative85_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio85_lower={log10_fraction(ratio85):.2f}")
    print(f"log10_best_cell_ratio85={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=85 rectangular certificate failed")


if __name__ == "__main__":
    main()
