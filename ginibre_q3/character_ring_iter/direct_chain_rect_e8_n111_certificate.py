#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=111.

This replays the n=113 rectangular machinery with N=111.  The positive
region is a 1/20 by 1/20 rectangular union inside |alpha.v| <= 4, and the
negative-region estimate is the same exact low-moment bound
176 * 16^(n-2).
"""

from __future__ import annotations

from fractions import Fraction

import direct_chain_rect_e8_n113_certificate as base
from direct_chain_a0_dominance import log10_fraction


N = 111
S_END = Fraction(537, 10)
EXPECTED_CELLS = 84702


def configure() -> None:
    base.N = N
    base.S_START = Fraction(40)
    base.S_END = S_END
    base.T_START = Fraction(30)
    base.T_END = Fraction(48)


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
    neg = base.negative_bound()

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

    ratio111 = sum(positive_terms) / neg
    best_term = max(positive_terms) / neg

    checks = {
        **cubic_checks,
        "sine_quartic_sturm": base.sine_scalar_check(),
        "root_identity_constants": base.root_identity_check(),
        "e_upper_1457_536": base.e_upper_check(),
        "low_moment_sources": base.low_moment_source_check(),
        "low_moment_negative_176": base.low_moment_negative_check(),
        f"cell_count_{EXPECTED_CELLS}": len(cells) == EXPECTED_CELLS,
        "all_cells_in_cap": all(
            s1 < base.cap_value() and t1 < base.cap_value() for _, s1, _, t1 in cells
        ),
        "positive_gaps": min(gaps) > 0,
        "positive_h111": min(h_values) > 0,
        "positive_direct111": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            base.exp_fractional_septic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            base.exp_fractional_septic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio111_gt_1": ratio111 > 1,
    }

    print(f"group={base.GROUP} n={N} delta={base.DELTA} e_upper={base.E_UPPER}")
    print(
        f"grid_step={base.GRID_STEP} s_start={base.S_START} s_end={base.S_END} "
        f"t_start={base.T_START} t_end={base.T_END}"
    )
    print(f"retained_cells={len(cells)} cap={base.cap_value()}")
    print(f"negative_bound={base.NEGATIVE_MOMENT_CONSTANT}*16^{N - 2}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap111_min={min(gaps)}")
    print(f"h111_min={min(h_values)} direct_factor111_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio111_lower={log10_fraction(ratio111):.2f}")
    print(f"log10_best_cell_ratio111={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=111 rectangular certificate failed")


if __name__ == "__main__":
    main()
