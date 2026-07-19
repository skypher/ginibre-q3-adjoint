#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=123.

This checks the arithmetic part of the DCT-RectTrig(E8, n=123) leaf.  The
region is a 1/40 by 1/40 rectangular union inside the root angle cap
|alpha.v| <= 4.  The local input is the same E8 quartic-sextic package,
with the gap coefficient 1/480 and sine denominator 20.

For the radial Gaussian and Weyl sine density losses, the checker uses the
fractional-part exponential lower bound from the n=129 bridge certificate:

    exp(-x) >= (1457/536)^(-floor(x))
              * (1 - r + r^2/2 - r^3/6),  r = x - floor(x).

The script compares the resulting rectangular-union lower bound for
D_E8(123) = Q_3^E8(125) - 4 Q_3^E8(123) with the negative-region bound
520*16^123.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, log10_fraction
from direct_chain_rect_e8_n131_certificate import root_identity_check


GROUP = "E8"
N = 123
D = 248
C_NEG = 8
RANK = 8
ALPHA = 250
SHAPE = 124
QUARTIC_DENOM = 50
SEXTIC_DENOM = 1800

DELTA = Fraction(4)
GRID_STEP = Fraction(1, 40)
S_START = Fraction(53)
S_END = Fraction(59)
T_START = Fraction(69, 2)
T_END = Fraction(75, 2)

PI_UPPER = Fraction(355, 113)
E_UPPER = Fraction(1457, 536)
GAP6_LOWER = Fraction(1, 480)
COS4_COEFF = Fraction(1, 12)
COS6_LOSS = Fraction(1, 360)
SINE_EXP_DENOM = 20


def cap_value() -> Fraction:
    # E8 has kappa = 30, so S_delta(n) = delta^2 * 30n/(4d).
    return DELTA * DELTA * 30 * N / (4 * D)


def e_upper_check() -> bool:
    # e < sum_{k=0}^7 1/k! + 1/(7!*7) < 1457/536.
    partial = sum(Fraction(1, factorial(k)) for k in range(8))
    return partial + Fraction(1, factorial(7) * 7) < E_UPPER


def exp_fractional_cubic_lower(value: Fraction) -> Fraction:
    whole = value.numerator // value.denominator
    frac = value - whole
    cubic = Fraction(1) - frac + frac * frac / 2 - frac**3 / 6
    return cubic / (E_UPPER**whole)


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    rho = Fraction(D, 6 * QUARTIC_DENOM * N)
    sixth_gain = GAP6_LOWER * Fraction((2 * D) ** 2, SEXTIC_DENOM)
    return s0 - t1 - rho * s1 * s1 + sixth_gain * s0**3 / (N * N)


def local_h_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> Fraction:
    quartic_gain = COS4_COEFF * Fraction(2 * D, QUARTIC_DENOM)
    sixth_loss = COS6_LOSS * Fraction((2 * D) ** 2, SEXTIC_DENOM)
    return (
        Fraction(1)
        - (s1 + t1) / N
        + quartic_gain * (s0 * s0 + t0 * t0) / (N * N)
        - sixth_loss * (s1**3 + t1**3) / (N**3)
    )


def local_direct_factor_lower(h: Fraction) -> Fraction:
    return Fraction((2 * D) ** 2) * h * h - 4


def candidate_cells() -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    cells: list[tuple[Fraction, Fraction, Fraction, Fraction]] = []
    cap = cap_value()
    s0 = S_START
    while s0 + GRID_STEP <= cap and s0 + GRID_STEP <= S_END:
        t0 = T_START
        while t0 + GRID_STEP <= T_END:
            s1 = s0 + GRID_STEP
            t1 = t0 + GRID_STEP
            if gap_lower(s0, s1, t0, t1) > 0:
                h = local_h_lower(s0, s1, t0, t1)
                if h > 0 and local_direct_factor_lower(h) > 0:
                    cells.append((s0, s1, t0, t1))
            t0 += GRID_STEP
        s0 += GRID_STEP
    return cells


def gamma_inv_square() -> Fraction:
    # Gamma(124) = 123!.
    return Fraction(1, factorial(SHAPE - 1) ** 2)


def rectangle_fraction_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction
) -> Fraction:
    s_integral = (s1**SHAPE - s0**SHAPE) / SHAPE
    t_integral = (t1**SHAPE - t0**SHAPE) / SHAPE
    return (
        gap_lower(s0, s1, t0, t1) ** 2
        * exp_fractional_cubic_lower(s1 + t1)
        * s_integral
        * t_integral
        * gamma_inv_square()
        / SHAPE
    )


def a0_lower() -> Fraction:
    return a0_prefactor(ROOT_DATA[GROUP]) / (PI_UPPER**RANK)


def sine_density_lower(s1: Fraction, t1: Fraction) -> Fraction:
    exponent = Fraction(4 * D, SINE_EXP_DENOM * N) * (s1 + t1)
    return exp_fractional_cubic_lower(exponent)


def negative_bound() -> int:
    return 2 * ((2 * C_NEG) ** 2 + 4) * (2 * C_NEG) ** N


def positive_rectangle_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, a0: Fraction
) -> Fraction:
    h = local_h_lower(s0, s1, t0, t1)
    direct = local_direct_factor_lower(h)
    if h <= 0 or direct <= 0:
        return Fraction(0)
    return (
        a0
        * Fraction((2 * D) ** N, N**ALPHA)
        * rectangle_fraction_lower(s0, s1, t0, t1)
        * h**N
        * direct
        * sine_density_lower(s1, t1)
    )


def ratio_lower(
    cells: list[tuple[Fraction, Fraction, Fraction, Fraction]], a0: Fraction
) -> Fraction:
    return sum(positive_rectangle_lower(*cell, a0) for cell in cells) / negative_bound()


def main() -> None:
    cells = candidate_cells()
    if len(cells) != 28121:
        raise AssertionError(f"expected 28121 retained grid cells, got {len(cells)}")
    gaps = [gap_lower(s0, s1, t0, t1) for s0, s1, t0, t1 in cells]
    h_values = [local_h_lower(s0, s1, t0, t1) for s0, s1, t0, t1 in cells]
    direct_values = [local_direct_factor_lower(h) for h in h_values]
    rect_fracs = [rectangle_fraction_lower(*cell) for cell in cells]
    sine_exponents = [
        Fraction(4 * D, SINE_EXP_DENOM * N) * (s1 + t1)
        for _, s1, _, t1 in cells
    ]
    a0 = a0_lower()
    ratio123 = ratio_lower(cells, a0)
    best_term = max(positive_rectangle_lower(*cell, a0) / negative_bound() for cell in cells)

    checks = {
        "root_identity_constants": root_identity_check(),
        "e_upper_1457_536": e_upper_check(),
        "cell_count_28121": len(cells) == 28121,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h123": min(h_values) > 0,
        "positive_direct123": min(direct_values) > 0,
        "positive_radial_exp_lower": all(
            exp_fractional_cubic_lower(s1 + t1) > 0 for _, s1, _, t1 in cells
        ),
        "positive_sine_exp_lower": all(
            exp_fractional_cubic_lower(exponent) > 0 for exponent in sine_exponents
        ),
        "ratio123_gt_1": ratio123 > 1,
    }

    print(f"group={GROUP} n={N} delta={DELTA} e_upper={E_UPPER}")
    print(
        f"grid_step={GRID_STEP} s_start={S_START} s_end={S_END} "
        f"t_start={T_START} t_end={T_END}"
    )
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap123_min={min(gaps)}")
    print(f"h123_min={min(h_values)} direct_factor123_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio123_lower={log10_fraction(ratio123):.2f}")
    print(f"log10_best_cell_ratio123={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=123 rectangular certificate failed")


if __name__ == "__main__":
    main()
