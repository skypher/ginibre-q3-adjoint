#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=107.

This replays the n=113 rectangular machinery with N=107.  The positive
region is a 1/20 by 1/20 rectangular union inside |alpha.v| <= 4.  The
negative-region estimate uses the asymmetric degree-28 moment majorant

    (63/65) * 16^(n-12) * S^12 * (S^2+4) * ((496-S)/512)^14

for S = X+Y on the negative-sign intervals S in [-16,-2] and S in [0,2].
"""

from __future__ import annotations

import re
from fractions import Fraction
from math import comb
from pathlib import Path

import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction
from direct_chain_rect_e8_n117_certificate import (
    poly_divmod,
    poly_eval,
    positive_no_roots,
)


N = 107
P_NEG = 12
ASYM_POWER = 14
S_END = Fraction(207, 4)
EXPECTED_CELLS = 70562
ASYMMETRIC_INTEGRAL = Fraction(
    9594616716049666093860048619540906615318349630451,
    21267647932558653966460912964485513216,
)
MOMENTS_0_30 = (
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
)


def configure() -> None:
    base.N = N
    base.S_START = Fraction(40)
    base.S_END = S_END
    base.T_START = Fraction(30)
    base.T_END = Fraction(48)
    base.GRID_STEP = Fraction(1, 20)


def parse_moment_prefix(path: Path) -> tuple[int, ...]:
    moments: dict[int, int] = {}
    with path.open() as handle:
        for line in handle:
            match = re.search(r"(?:^|\bm_)\s*(\d+)\s*(?:=)?\s+(-?\d+)", line)
            if match:
                moments[int(match.group(1))] = int(match.group(2))
    return tuple(moments[k] for k in range(len(MOMENTS_0_30)))


def moment_source_check() -> bool:
    here = Path(__file__).resolve().parent
    expected = tuple(MOMENTS_0_30)
    return (
        parse_moment_prefix(here / "logs" / "e8_70.log") == expected
        and parse_moment_prefix(here.parent / "references" / "oeis_A179663.txt") == expected
    )


def moment_integral(power: int) -> int:
    total = 0
    for k in range(power + 1):
        coeff = comb(power, k)
        total += coeff * MOMENTS_0_30[k + 2] * MOMENTS_0_30[power - k]
        total += coeff * MOMENTS_0_30[k] * MOMENTS_0_30[power - k + 2]
        total -= 2 * coeff * MOMENTS_0_30[k + 1] * MOMENTS_0_30[power - k + 1]
    return total


def poly_trim(poly: list[Fraction]) -> list[Fraction]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    result = [Fraction(0)] * size
    for i in range(size):
        result[i] = (left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0)
    return poly_trim(result)


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_coeff in enumerate(left):
        for j, right_coeff in enumerate(right):
            result[i + j] += left_coeff * right_coeff
    return poly_trim(result)


def poly_scale(poly: list[Fraction], factor: int | Fraction) -> list[Fraction]:
    return poly_trim([factor * coeff for coeff in poly])


def poly_pow(poly: list[Fraction], exponent: int) -> list[Fraction]:
    result = [Fraction(1)]
    for _ in range(exponent):
        result = poly_mul(result, poly)
    return result


def asymmetric_negative_interval_check() -> bool:
    # With x=-S, check the residual after clearing 65*512^14 and then
    # dividing by the endpoint factor 16-x.
    x = [Fraction(0), Fraction(1)]
    x2 = poly_mul(x, x)
    lhs = poly_scale(
        poly_mul(poly_add(x2, [Fraction(4)]), poly_pow([Fraction(496), Fraction(1)], 14)),
        63 * 16 ** (N - P_NEG),
    )
    rhs = poly_scale(
        poly_mul(poly_pow(x, N - P_NEG), poly_add(x2, [Fraction(-4)])),
        65 * 512**14,
    )
    residual = poly_add(lhs, poly_scale(rhs, -1))
    quotient, remainder = poly_divmod(residual, [Fraction(16), Fraction(-1)])
    return (
        not remainder
        and poly_eval(quotient, Fraction(2)) > 0
        and poly_eval(quotient, Fraction(16)) > 0
        and positive_no_roots(quotient, Fraction(2), Fraction(16))
    )


def asymmetric_positive_interval_check() -> bool:
    # On 0 <= S <= 2, divide the residual by S^12; the omitted factor is
    # non-negative and both sides vanish at S=0.
    s = [Fraction(0), Fraction(1)]
    s2 = poly_mul(s, s)
    lhs = poly_scale(
        poly_mul(poly_add(s2, [Fraction(4)]), poly_pow([Fraction(496), Fraction(-1)], 14)),
        63 * 16 ** (N - P_NEG),
    )
    rhs = poly_scale(
        poly_mul(poly_pow(s, N - P_NEG), poly_add([Fraction(4)], poly_scale(s2, -1))),
        65 * 512**14,
    )
    residual = poly_add(lhs, poly_scale(rhs, -1))
    return positive_no_roots(residual, Fraction(0), Fraction(2))


def asymmetric_integral() -> Fraction:
    total = Fraction(0)
    for i in range(ASYM_POWER + 1):
        coeff = comb(ASYM_POWER, i) * 496 ** (ASYM_POWER - i) * (-1) ** i
        total += Fraction(coeff, 512**ASYM_POWER) * (
            moment_integral(P_NEG + 2 + i) + 4 * moment_integral(P_NEG + i)
        )
    return total


def asymmetric_moment_check() -> bool:
    return (
        moment_integral(12) == 5428655636
        and moment_integral(14) == 1070740575644
        and asymmetric_integral() == ASYMMETRIC_INTEGRAL
    )


def negative_bound() -> Fraction:
    return Fraction(63, 65) * 16 ** (N - P_NEG) * ASYMMETRIC_INTEGRAL


def old_negative_bound() -> int:
    return (moment_integral(14) + 4 * moment_integral(12)) * 16 ** (N - P_NEG)


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

    ratio107 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg

    checks = {
        **cubic_checks,
        "sine_quartic_sturm": base.sine_scalar_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        "moment_sources_0_30": moment_source_check(),
        "asymmetric_negative_interval_sturm": asymmetric_negative_interval_check(),
        "asymmetric_positive_interval_sturm": asymmetric_positive_interval_check(),
        "asymmetric_moment_integral": asymmetric_moment_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(
            s1 < base.cap_value() and t1 < base.cap_value() for _, s1, _, t1 in cells
        ),
        "positive_gaps": min(gaps) > 0,
        "positive_h107": min(h_values) > 0,
        "positive_direct107": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio107_gt_1": ratio107 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={base.DELTA} e_upper={base.E_UPPER}")
    print(
        f"grid_step={base.GRID_STEP} s_start={base.S_START} s_end={base.S_END} "
        f"t_start={base.T_START} t_end={base.T_END}"
    )
    print(f"retained_cells={len(cells)} cap={base.cap_value()}")
    print(f"asym_integral={ASYMMETRIC_INTEGRAL}")
    print(f"log10_asym_negative_vs_p12={log10_fraction(neg / old_negative_bound()):.2f}")
    print("negative_bound=(63/65)*16^95*asym_integral")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap107_min={min(gaps)}")
    print(f"h107_min={min(h_values)} direct_factor107_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio107_lower={log10_fraction(ratio107):.2f}")
    print(f"log10_best_cell_ratio107={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=107 rectangular certificate failed")


if __name__ == "__main__":
    main()
