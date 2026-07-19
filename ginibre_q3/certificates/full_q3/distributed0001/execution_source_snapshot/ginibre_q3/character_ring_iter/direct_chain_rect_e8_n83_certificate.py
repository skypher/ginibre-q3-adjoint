#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=83.

This reuses the n=85 rectangular machinery with N=83.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_92, checked against the local m_0..m_70 character-ring log.
This permits a degree-38 square majorant at the same p=12 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n89_certificate as template
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n99_n97_certificate as shared
from direct_chain_a0_dominance import log10_fraction


N = 83
P_NEG = 12
NEG_DEGREE = 38
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 91138
MOMENTS_0_92 = prev95.MOMENTS_0_70 + (
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
    7566385447785130244097480106783651234709147512629865928330665298909382359168224089309482199841301279446987381885863040618855936,
    643539682611983635221052667467181610004835846072031095557465137562914783242754873994458040809944307301737283417986291313825250384,
    55243536790553384757856474131522285658424037485294992240244823864378023908294530262574782229751671326703212065765182837042785328912,
    4785608534006322970629517207008860446800991091514759387399534753159452436410778964305557322373041621673115866259169663944528197154888,
    418287144523196219371866155170789986707226492225704076967013598709841498794913059918976450530726410902295409765563712045308089767188896,
    36883096437939070423637216037426824097540846213687258071316706754324839822776931784481118499591302325505163049736816327363461025875969392,
    3280433922779716659235757169218135489448361965928393342322481519338804120878949142603599891244921337737510899532169778813444559330767388504,
    294254741448395650521947765533072609866216719852588539783798087814740641036511546090999119522229986774987753267396048498283526176339524056284,
)

NEG_CHEB_COEFFS = (
    Fraction("19442497235669361/62500000000000000000"),
    Fraction("-266297900053345251/1000000000000000000000"),
    Fraction("122565762229486143/1000000000000000000000"),
    Fraction("4670020817725957/31250000000000000000"),
    Fraction("-288135915501630741/500000000000000000000"),
    Fraction("444893053682977413/500000000000000000000"),
    Fraction("-50172369622834359/625000000000000000000"),
    Fraction("-339695175467239689/500000000000000000000"),
    Fraction("161390172531801371/250000000000000000000"),
    Fraction("-355340495291895709/1000000000000000000000"),
    Fraction("35466689688532501/2500000000000000000000"),
    Fraction("295736673508335747/1000000000000000000000"),
    Fraction("-98961611190745923/200000000000000000000"),
    Fraction("9785381029456463/20000000000000000000"),
    Fraction("-216087515267137991/1000000000000000000000"),
    Fraction("-249886904213621837/1000000000000000000000"),
    Fraction("630471234658398793/1000000000000000000000"),
    Fraction("-309134570488522713/500000000000000000000"),
    Fraction("359125183600773/2000000000000000000"),
    Fraction("24679459118226217/62500000000000000000"),
    Fraction("-363063741084753951/500000000000000000000"),
    Fraction("293476722237992211/500000000000000000000"),
    Fraction("-69373390750091573/2000000000000000000000"),
    Fraction("-117068597276593733/200000000000000000000"),
    Fraction("819763401645266787/1000000000000000000000"),
    Fraction("-209841318096667197/500000000000000000000"),
    Fraction("-182664516282605397/500000000000000000000"),
    Fraction("43932168397696757/50000000000000000000"),
    Fraction("-603062632060651343/1000000000000000000000"),
    Fraction("-260101134681331621/1000000000000000000000"),
    Fraction("833401858016795849/1000000000000000000000"),
    Fraction("-500456730427659689/1000000000000000000000"),
    Fraction("-65385674507054993/200000000000000000000"),
    Fraction("667803915578651973/1000000000000000000000"),
    Fraction("-41473305706972731/250000000000000000000"),
    Fraction("-22846251637431727/50000000000000000000"),
    Fraction("31456980636188547/100000000000000000000"),
    Fraction("13583722755901247/31250000000000000000"),
    Fraction("56868794477882047/500000000000000000000"),
)

NEGATIVE_INTEGRAL = Fraction(
    69506923908375742396638750904720698654253961502478284262296998242926162620300793054819096177856546091509678268326171062873371758355860924644968831631635906940012075850490819786026964402425828891779371420933273301225535316551879001356040185988627,
    61011724716641032694863951249238379346746835550611863791471897122554948049358208648682622736448229569001250560659210985228745906718984537731398722994440824073270393246860914952526114334982958037512832444025000000000000000000000000000000000000000000,
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
    template.MOMENTS_0_72 = MOMENTS_0_92
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
        name.replace("n89", "n83"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    arxiv_moments: dict[int, int] = {}
    with (here.parent / "references" / "arxiv_2412_21189_e8_m71_m92.txt").open() as handle:
        for line in handle:
            match = shared.re.search(r"\bm_\s*(\d+)\s*=\s+(-?\d+)", line)
            if match:
                arxiv_moments[int(match.group(1))] = int(match.group(2))
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
        and arxiv_moments == {k: MOMENTS_0_92[k] for k in range(71, 93)}
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

    ratio83 = positive_total / neg
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
        "positive_h83": min(h_values) > 0,
        "positive_direct83": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio83_gt_1": ratio83 > 1,
    }

    print(f"group={template.base.GROUP} n={N} delta={DELTA} e_upper={template.base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^71*optimized_degree38_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap83_min={min(gaps)}")
    print(f"h83_min={min(h_values)} direct_factor83_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive83_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative83_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio83_lower={log10_fraction(ratio83):.2f}")
    print(f"log10_best_cell_ratio83={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=83 rectangular certificate failed")


if __name__ == "__main__":
    main()
