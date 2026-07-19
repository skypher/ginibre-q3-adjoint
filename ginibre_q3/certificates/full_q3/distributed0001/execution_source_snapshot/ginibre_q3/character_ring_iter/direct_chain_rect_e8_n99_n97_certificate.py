#!/usr/bin/env python3
"""Exact-rational certificate for the E8 bridge steps at n=99 and n=97.

This replays the n=101 rectangular machinery with one shared shifted
Chebyshev degree-26 negative-region majorant

    (63/65) * 16^(n-8) * S^8 * (S^2+4) * C_26(S)^2,
    C_26(S) = T_26((S-249)/247) / T_26(-265/247),

for S = X+Y on the negative-sign intervals S in [-16,-2] and S in [0,2].
The interval residuals are checked by exact Bernstein-coefficient positivity
after clearing the non-negative endpoint factors.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from math import comb
from pathlib import Path

import direct_chain_rect_e8_n105_certificate as prev
import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import poly_divmod, poly_eval


@dataclass(frozen=True)
class BridgeCase:
    n: int
    s_end: Fraction
    t_end: Fraction
    expected_cells: int


P_NEG = 8
CHEB_DEGREE = 26
MOMENT_MAX = P_NEG + 2 * CHEB_DEGREE + 4
BRIDGE_CASES = (
    BridgeCase(n=99, s_end=Fraction(479, 10), t_end=Fraction(479, 10), expected_cells=42662),
    BridgeCase(n=97, s_end=Fraction(469, 10), t_end=Fraction(469, 10), expected_cells=35869),
)
CHEBYSHEV_INTEGRAL = Fraction(
    652708718049108501279891206610792628257515262318580914182134886342747101135881782357617401218658855060065796860047212614264650354628,
    2445200692443756455056633566892545977519407487004403208271763856863104974462205280301950124320130386149081237313785855952197768855841,
)
MOMENTS_0_64 = (
    1,
    0,
    1,
    1,
    5,
    16,
    79,
    421,
    2674,
    19244,
    156612,
    1423028,
    14320350,
    158390872,
    1912977222,
    25083283995,
    355246037162,
    5409471180024,
    88200546561838,
    1534120589972637,
    28369229081383675,
    556021169447494656,
    11517512836906556032,
    251487262264563372960,
    5774531032488992734964,
    139120785179133253295792,
    3509477193620539275437180,
    92518968639950402918894032,
    2544345801307983812697346302,
    72869933800796610152638668712,
    2170019149477016767674782447224,
    67093405561880047721756198386919,
    2150768850578662380494169570944344,
    71389941002645611223406888140851312,
    2450607577121624281936293229672316980,
    86895607835705190487787330116296398491,
    3179292476229902262666013065673766456393,
    119899791672381439086884351990618657158992,
    4656202599065199614841650099247510706257419,
    186021369312670669651530754293447277895990825,
    7638774001497058254021480118216385880012703810,
    322140681015850959972608029794754059375800591564,
    13940470504581099713887149778641872549235207642284,
    618563901983246708306382342380114037201295177009956,
    28122099078999983794691251963812813407439800694359474,
    1309067661627492536592017480200904551886019921861594056,
    62349987400759433187192618148879612311291607505176226994,
    3036630811208529343942589627171060237916769491750479631533,
    151134063097209505953250239274102960090325015819074565743106,
    7682314964369183858898705928934705343540905116367973435660888,
    398599271660325035613108964432283517489951642195459285940349252,
    21098947726823578520650375872360786720699889109330152834107533330,
    1138779671661759560160327745985285364682157893025942844565601725142,
    62640739094612404528291453697528335794193709593483751013645630950432,
    3509977977832368731937441900497516562476975797535341264129377207968146,
    200254533368267627610283450566048163803321402532169437126989985226578638,
    11627823930341636925131668617520986536644617221659519936125267797819665876,
    686858909943128658185391537823511101669776180401218887150970949832957611128,
    41258299679392448426937932545773814197475129151374453127809553948400286903340,
    2519173073447130141383736879933095769286172232232795384946790480383259571217584,
    156293865552003246895617231956808617360541851719355695877503247080326092883833200,
    9849253890124530737954155151827574115900678109277195872991241895991661053656369792,
    630213960291872584200180205168892664212008403525291694395502597987923480090222752570,
    40930507452224784448386263411352077381363161137223130745566919000735207420042669316928,
    2697347333459503776603530486147103258522010653168049200795738727849810852537162110102072,
)


def configure(case: BridgeCase) -> None:
    base.N = case.n
    base.S_START = Fraction(40)
    base.S_END = case.s_end
    base.T_START = Fraction(30)
    base.T_END = case.t_end
    base.GRID_STEP = Fraction(1, 20)


def parse_moment_prefix(path: Path, length: int) -> tuple[int, ...]:
    moments: dict[int, int] = {}
    with path.open() as handle:
        for line in handle:
            match = re.search(r"(?:^|\bm_)\s*(\d+)\s*(?:=)?\s+(-?\d+)", line)
            if match:
                moments[int(match.group(1))] = int(match.group(2))
    return tuple(moments[k] for k in range(length))


def moment_source_check() -> bool:
    here = Path(__file__).resolve().parent
    return (
        parse_moment_prefix(here / "logs" / "e8_70.log", len(MOMENTS_0_64))
        == MOMENTS_0_64
        and parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt", 31)
        == MOMENTS_0_64[:31]
    )


def moment_integral(power: int) -> int:
    total = 0
    for k in range(power + 1):
        coeff = comb(power, k)
        total += coeff * MOMENTS_0_64[k + 2] * MOMENTS_0_64[power - k]
        total += coeff * MOMENTS_0_64[k] * MOMENTS_0_64[power - k + 2]
        total -= 2 * coeff * MOMENTS_0_64[k + 1] * MOMENTS_0_64[power - k + 1]
    return total


def chebyshev_factor() -> list[Fraction]:
    raw = prev.compose_linear(
        prev.chebyshev_polynomial(CHEB_DEGREE),
        Fraction(1, 247),
        Fraction(-249, 247),
    )
    value_at_minus_16 = poly_eval(raw, Fraction(-16))
    return [coeff / value_at_minus_16 for coeff in raw]


def compose_interval(poly: list[Fraction], left: Fraction, right: Fraction) -> list[Fraction]:
    return prev.compose_linear(poly, right - left, left)


def bernstein_coefficients_on_unit_interval(poly: list[Fraction]) -> list[Fraction]:
    degree = len(poly) - 1
    coefficients: list[Fraction] = []
    for i in range(degree + 1):
        value = Fraction(0)
        for k in range(i + 1):
            if k < len(poly):
                value += poly[k] * Fraction(comb(i, k), comb(degree, k))
        coefficients.append(value)
    return coefficients


def bernstein_positive_on(poly: list[Fraction], left: Fraction, right: Fraction) -> bool:
    unit_poly = compose_interval(poly, left, right)
    return all(coeff > 0 for coeff in bernstein_coefficients_on_unit_interval(unit_poly))


def chebyshev_negative_interval_check(n: int) -> bool:
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
        Fraction(63, 65) * 16 ** (n - P_NEG),
    )
    rhs = prev.helper.poly_mul(prev.helper.poly_pow(x, n), prev.helper.poly_add(x2, [Fraction(-4)]))
    quotient, remainder = poly_divmod(prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1)), xp)
    quotient2, remainder2 = poly_divmod(quotient, [Fraction(16), Fraction(-1)])
    return (
        not remainder
        and not remainder2
        and poly_eval(quotient2, Fraction(2)) > 0
        and poly_eval(quotient2, Fraction(16)) > 0
        and bernstein_positive_on(quotient2, Fraction(2), Fraction(16))
    )


def chebyshev_positive_interval_check(n: int) -> bool:
    s_poly = chebyshev_factor()
    s = [Fraction(0), Fraction(1)]
    s2 = prev.helper.poly_mul(s, s)
    sp = prev.helper.poly_pow(s, P_NEG)
    lhs = prev.helper.poly_scale(
        prev.helper.poly_mul(
            prev.helper.poly_mul(sp, prev.helper.poly_add(s2, [Fraction(4)])),
            prev.helper.poly_mul(s_poly, s_poly),
        ),
        Fraction(63, 65) * 16 ** (n - P_NEG),
    )
    rhs = prev.helper.poly_mul(
        prev.helper.poly_pow(s, n),
        prev.helper.poly_add([Fraction(4)], prev.helper.poly_scale(s2, -1)),
    )
    quotient, remainder = poly_divmod(prev.helper.poly_add(lhs, prev.helper.poly_scale(rhs, -1)), sp)
    return (
        not remainder
        and poly_eval(quotient, Fraction(0)) > 0
        and poly_eval(quotient, Fraction(2)) > 0
        and bernstein_positive_on(quotient, Fraction(0), Fraction(2))
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


def negative_bound(n: int) -> Fraction:
    return Fraction(63, 65) * 16 ** (n - P_NEG) * CHEBYSHEV_INTEGRAL


def old_negative_bound(n: int) -> int:
    return (moment_integral(14) + 4 * moment_integral(12)) * 16 ** (n - 12)


def case_ratios(case: BridgeCase) -> tuple[dict[str, bool], dict[str, Fraction | int]]:
    configure(case)
    cells = base.candidate_cells()
    if len(cells) != case.expected_cells:
        raise AssertionError(
            f"n={case.n}: expected {case.expected_cells} retained grid cells, got {len(cells)}"
        )

    gaps: list[Fraction] = []
    h_values: list[Fraction] = []
    direct_values: list[Fraction] = []
    rect_fracs: list[Fraction] = []
    positive_terms: list[Fraction] = []
    sine_exponents: list[Fraction] = []
    a0 = base.a0_lower()
    neg = negative_bound(case.n)

    for cell in cells:
        s0, s1, t0, t1 = cell
        gap = base.gap_lower(s0, s1, t0, t1)
        h = base.local_h_lower(s0, s1, t0, t1)
        direct = base.local_direct_factor_lower(h)
        rect_frac = base.rectangle_fraction_lower(s0, s1, t0, t1)
        sine_exponent = (
            2 * base.SINE_QUADRATIC * Fraction(2 * base.D) * (s1 + t1) / case.n
            + 2
            * base.SINE_QUARTIC
            * Fraction((2 * base.D) ** 2, base.QUARTIC_DENOM)
            * (s1 * s1 + t1 * t1)
            / (case.n * case.n)
        )
        sine = base.exp_fractional_septic_lower(sine_exponent)
        term = (
            a0
            * Fraction((2 * base.D) ** case.n, case.n**base.ALPHA)
            * rect_frac
            * h**case.n
            * direct
            * sine
        )
        gaps.append(gap)
        h_values.append(h)
        direct_values.append(direct)
        rect_fracs.append(rect_frac)
        positive_terms.append(term)
        sine_exponents.append(sine_exponent)

    ratio = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg
    checks = {
        f"cell_count_{case.n}_{case.expected_cells}": len(cells) == case.expected_cells,
        f"all_cells_in_cap_{case.n}": all(
            s1 < base.cap_value() and t1 < base.cap_value() for _, s1, _, t1 in cells
        ),
        f"positive_gaps_{case.n}": min(gaps) > 0,
        f"positive_h_{case.n}": min(h_values) > 0,
        f"positive_direct_{case.n}": min(direct_values) > 0,
        f"positive_radial_exp_lower_{case.n}": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        f"positive_sine_exp_lower_{case.n}": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        f"ratio{case.n}_gt_1": ratio > 1,
    }
    data: dict[str, Fraction | int] = {
        "cells": len(cells),
        "cap": base.cap_value(),
        "union_fraction": sum(rect_fracs),
        "gap_min": min(gaps),
        "h_min": min(h_values),
        "direct_min": min(direct_values),
        "sine_exponent_max": max(sine_exponents),
        "ratio": ratio,
        "best_term": best_term,
        "negative_vs_p12": neg / old_negative_bound(case.n),
    }
    return checks, data


def main() -> None:
    cubic_checks = base.scalar_polynomial_checks()
    cheb = chebyshev_factor()
    shared_checks = {
        **cubic_checks,
        "sine_quartic_sturm": base.sine_scalar_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        f"moment_log_0_{MOMENT_MAX}_oeis_overlap_0_30": moment_source_check(),
        "shifted_chebyshev_normalization": poly_eval(cheb, Fraction(-16)) == 1,
        "shifted_chebyshev_moment_integral": chebyshev_moment_check(),
    }

    all_checks: dict[str, bool] = dict(shared_checks)
    case_data: dict[int, dict[str, Fraction | int]] = {}
    for case in BRIDGE_CASES:
        all_checks[f"shifted_chebyshev_negative_bernstein_{case.n}"] = (
            chebyshev_negative_interval_check(case.n)
        )
        all_checks[f"shifted_chebyshev_positive_bernstein_{case.n}"] = (
            chebyshev_positive_interval_check(case.n)
        )
        checks, data = case_ratios(case)
        all_checks.update(checks)
        case_data[case.n] = data

    print(
        f"group={base.GROUP} n_values={','.join(str(case.n) for case in BRIDGE_CASES)} "
        f"delta={base.DELTA} e_upper={base.E_UPPER}"
    )
    print(
        f"grid_step={base.GRID_STEP} s_start=40 t_start=30 "
        f"shifted_chebyshev_degree={CHEB_DEGREE} moment_max={MOMENT_MAX}"
    )
    print(f"cheb_at_0={poly_eval(cheb, Fraction(0))}")
    print(f"cheb_integral={CHEBYSHEV_INTEGRAL}")
    for case in BRIDGE_CASES:
        data = case_data[case.n]
        print(f"case n={case.n} s_end={case.s_end} t_end={case.t_end}")
        print(f"retained_cells={data['cells']} cap={data['cap']}")
        print(f"negative_bound=(63/65)*16^{case.n - P_NEG}*cheb_integral")
        print(f"log10_cheb_negative_vs_p12={log10_fraction(data['negative_vs_p12']):.2f}")
        print(f"log10_union_fraction_lower={log10_fraction(data['union_fraction']):.2f}")
        print(f"gap{case.n}_min={data['gap_min']}")
        print(f"h{case.n}_min={data['h_min']} direct_factor{case.n}_min={data['direct_min']}")
        print(f"sine_exponent{case.n}_max={data['sine_exponent_max']}")
        print(f"log10_ratio{case.n}_lower={log10_fraction(data['ratio']):.2f}")
        print(f"log10_best_cell_ratio{case.n}={log10_fraction(data['best_term']):.2f}")
    for name, ok in all_checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(all_checks.values()):
        raise AssertionError("E8 n=99,n=97 rectangular range certificate failed")


if __name__ == "__main__":
    main()
