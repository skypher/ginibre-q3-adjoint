#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=81.

This reuses the n=83 rectangular machinery with N=81.  The positive side
keeps the root-angle cap |alpha.v| <= 14/3, grid step 1/20, the cap-adjusted
cubic gap coefficient, and the sextic Weyl-sine scalar bound.

The negative side uses the arXiv 2412.21189 ancillary E8 adjoint tensor-rank
values m_71,...,m_100, checked against the local m_0..m_70 character-ring log.
This permits a degree-42 square majorant at the same p=12 moment weight.
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n89_certificate as template
import direct_chain_rect_e8_n95_certificate as prev95
import direct_chain_rect_e8_n99_n97_certificate as shared
from direct_chain_a0_dominance import log10_fraction


N = 81
P_NEG = 12
NEG_DEGREE = 42
MOMENT_MAX = P_NEG + 2 * NEG_DEGREE + 4
DELTA = Fraction(14, 3)
GRID_STEP = Fraction(1, 20)
S_START = Fraction(40)
T_START = Fraction(30)
EXPECTED_CELLS = 79060
NEGATIVE_BERNSTEIN_PIECES = 8
MOMENTS_0_100 = prev95.MOMENTS_0_70 + (
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
    26615993029156185706819187624210902117084829191506053914949138689287637800017519244852663031430153214445558816120219787540939990890610469401328,
    2427332166680152352118369722655067011113158551489731266591855888524387501935533975095484204431393480034463370819074768221080161071528701155209628,
    223164372309640260500083436346711257062400763522820657142937080083654288692033007655228674601826886343916413552108663764041827935593399713577352170,
    20681061956354021357368389680134284231581284236875918573927975301645168618130176600974903921139326754317987689411777977451689235819480079172231780000,
    1931602350600318284893371178345070735477182448529952712443931522480426332617029689552991991046942009346541114427350036091265483358729751826425059864928,
    181804847792418636498167299450240940363227051339726242195555386759273650183719576073998739703113834220124660905422193043733256525286129744961191596924136,
    17241826130158021565221742880478393608102624030314803878163596056014570703401234755487523496761414724809547895624581604695267499905790032317347139181650147,
    1647402066054989633758123510269365166279880064500442045060407411700348494535823374334785155244179858975641874327973274326670123750796213931716347633680512809,
)

NEG_CHEB_COEFFS = (
    Fraction("-150157430179965833/1000000000000000000000"),
    Fraction("63977911644115559/500000000000000000000"),
    Fraction("-284942843027740671/5000000000000000000000"),
    Fraction("-47589662612630537/625000000000000000000"),
    Fraction("281768354692001847/1000000000000000000000"),
    Fraction("-424075490137007581/1000000000000000000000"),
    Fraction("112521741825564263/5000000000000000000000"),
    Fraction("334069079389622903/1000000000000000000000"),
    Fraction("-306213756806484889/1000000000000000000000"),
    Fraction("160204098300831071/1000000000000000000000"),
    Fraction("243084391603139077/50000000000000000000000"),
    Fraction("-74757888515537149/500000000000000000000"),
    Fraction("47295985521653383/200000000000000000000"),
    Fraction("-44780305900792483/200000000000000000000"),
    Fraction("1753725455950569/20000000000000000000"),
    Fraction("32528833706065763/250000000000000000000"),
    Fraction("-148501236185712691/500000000000000000000"),
    Fraction("277825877131495627/1000000000000000000000"),
    Fraction("-348719328388132733/5000000000000000000000"),
    Fraction("-188798633292219179/1000000000000000000000"),
    Fraction("65941676477006847/200000000000000000000"),
    Fraction("-16276205107823461/62500000000000000000"),
    Fraction("163721612356304227/10000000000000000000000"),
    Fraction("31822351943847439/125000000000000000000"),
    Fraction("-359264722772502691/1000000000000000000000"),
    Fraction("101033575891996491/500000000000000000000"),
    Fraction("7584605686212761/62500000000000000000"),
    Fraction("-35981762079918513/100000000000000000000"),
    Fraction("18498967265721459/62500000000000000000"),
    Fraction("25382286251133153/1250000000000000000000"),
    Fraction("-75075443530495719/250000000000000000000"),
    Fraction("277391322262920221/1000000000000000000000"),
    Fraction("-115329515380060077/10000000000000000000000"),
    Fraction("-206848357379179563/1000000000000000000000"),
    Fraction("24214256708592489/125000000000000000000"),
    Fraction("-1029787944194589/78125000000000000000"),
    Fraction("-130564219576533963/1000000000000000000000"),
    Fraction("57047728242605987/500000000000000000000"),
    Fraction("52034629665554037/5000000000000000000000"),
    Fraction("-941162028758150157/10000000000000000000000"),
    Fraction("396718652080757303/10000000000000000000000"),
    Fraction("794194112266904831/10000000000000000000000"),
    Fraction("224380753364137949/10000000000000000000000"),
)

NEGATIVE_INTEGRAL = Fraction(
    162988295159444920217175613977736273000922799417640375650345408112354141693199993604285605069801764134398766069421665451134117097220188593300942224975121288698309839236172476378298497541672298842570954921593146535995554007868966392982351229089785126310293950485983753,
    1662275213223649722955032911709714874544854221093257063513796535114453366206250958407904534456796140022503212421180405346282041850376672598241125169818513010264936158883955387393292011765444675330836543152256442397205420466966875000000000000000000000000000000000000000000,
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
    template.MOMENTS_0_72 = MOMENTS_0_100
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
        name.replace("n89", "n81"): ok
        for name, ok in template.scalar_polynomial_checks().items()
    }


def moment_source_check() -> bool:
    here = shared.Path(__file__).resolve().parent
    arxiv_moments: dict[int, int] = {}
    with (here.parent / "references" / "arxiv_2412_21189_e8_m71_m100.txt").open() as handle:
        for line in handle:
            match = shared.re.search(r"\bm_\s*(\d+)\s*=\s+(-?\d+)", line)
            if match:
                arxiv_moments[int(match.group(1))] = int(match.group(2))
    return (
        shared.parse_moment_prefix(here / "logs" / "e8_70.log", len(prev95.MOMENTS_0_70))
        == prev95.MOMENTS_0_70
        and shared.parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == prev95.MOMENTS_0_70[:31]
        and arxiv_moments == {k: MOMENTS_0_100[k] for k in range(71, 101)}
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


def optimized_negative_subdivided_interval_check() -> bool:
    install_parameters()
    factor = template.negative_factor()
    x = [Fraction(0), Fraction(1)]
    x2 = template.prev.helper.poly_mul(x, x)
    p_neg_x = template.prev.compose_linear(factor, Fraction(-1), Fraction(0))
    lhs = template.prev.helper.poly_scale(
        template.prev.helper.poly_mul(
            template.prev.helper.poly_add(x2, [Fraction(4)]),
            template.prev.helper.poly_mul(p_neg_x, p_neg_x),
        ),
        16 ** (N - P_NEG),
    )
    rhs = template.prev.helper.poly_mul(
        template.prev.helper.poly_pow(x, N - P_NEG),
        template.prev.helper.poly_add(x2, [Fraction(-4)]),
    )
    residual = template.prev.helper.poly_add(lhs, template.prev.helper.poly_scale(rhs, -1))
    width = Fraction(14, NEGATIVE_BERNSTEIN_PIECES)
    return (
        template.poly_eval(residual, Fraction(2)) > 0
        and template.poly_eval(residual, Fraction(16)) > 0
        and all(
            shared.bernstein_positive_on(
                residual, Fraction(2) + i * width, Fraction(2) + (i + 1) * width
            )
            for i in range(NEGATIVE_BERNSTEIN_PIECES)
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

    ratio81 = positive_total / neg
    checks = {
        **scalar_checks,
        "sine_sextic_14_3_sturm": template.sine_sextic_check(),
        "root_identity_constants": template.base.root_identity_check(),
        "e_upper_1457_536": template.base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_arxiv_extension": moment_source_check(),
        f"optimized_negative_interval_bernstein_{NEGATIVE_BERNSTEIN_PIECES}": (
            optimized_negative_subdivided_interval_check()
        ),
        "optimized_positive_interval_bernstein": template.optimized_positive_interval_check(),
        "optimized_negative_moment_integral": optimized_negative_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h81": min(h_values) > 0,
        "positive_direct81": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            template.base.exp_fractional_septic_lower(s1 + t1) > 0
            for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            template.base.exp_fractional_septic_lower(exponent) > 0
            for exponent in sine_exponents
        ),
        "ratio81_gt_1": ratio81 > 1,
    }

    print(f"group={template.base.GROUP} n={N} delta={DELTA} e_upper={template.base.E_UPPER}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print("gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000")
    print("sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)")
    print(f"negative_majorant_degree={NEG_DEGREE} shifted_basis_interval=[7,496]")
    print(f"negative_interval_subdivision_count={NEGATIVE_BERNSTEIN_PIECES}")
    print(f"negative_integral={negative_integral()}")
    print(f"log10_optimized_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=16^69*optimized_degree42_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap81_min={min(gaps)}")
    print(f"h81_min={min(h_values)} direct_factor81_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_positive81_lower={log10_fraction(positive_total):.2f}")
    print(f"log10_negative81_upper={log10_fraction(neg):.2f}")
    print(f"log10_ratio81_lower={log10_fraction(ratio81):.2f}")
    print(f"log10_best_cell_ratio81={log10_fraction(best_term / neg):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=81 rectangular certificate failed")


if __name__ == "__main__":
    main()
