#!/usr/bin/env python3
"""Exact-rational certificate for the E8 direct-Chain bridge step at n=131.

This checks the arithmetic part of the DCT-RectTrig(E8, n=131) leaf.  The
region is a 1/40 by 1/40 rectangular union inside the slightly wider root
angle cap |alpha.v| <= 19/6:

    s0 in {37, 1481/40, ..., floor-cap cells},
    t0 in {26, 1041/40, ..., 1239/40}.

The local input is the E8 quartic and sextic root identities

    sum_{alpha>0} (alpha.u)^4 = Q(u)^2 / 50,
    sum_{alpha>0} (alpha.u)^6 = Q(u)^3 / 1800.

The checker uses rational one-sided transcendental bounds

    pi < 355/113,     e < 11/4,

and the scalar inequalities on |z| <= 19/6

    2(1-cos z) >= z^2 - z^4/12 + z^6/432,
    2 cos z >= 2 - z^2 + z^4/12 - z^6/360,
    sin(z/2)/(z/2) >= exp(-z^2/21).

The script compares the resulting rectangular-union lower bound for
D_E8(131) = Q_3^E8(133) - 4 Q_3^E8(131) with the negative-region bound
520*16^131.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import factorial

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, log10_fraction


GROUP = "E8"
N = 131
D = 248
C_NEG = 8
RANK = 8
ALPHA = 250
SHAPE = 124
QUARTIC_DENOM = 50
SEXTIC_DENOM = 1800

DELTA = Fraction(19, 6)
GRID_STEP = Fraction(1, 40)
S_START = Fraction(37)
T_START = Fraction(26)
T_END = Fraction(31)

PI_UPPER = Fraction(355, 113)
E_UPPER = Fraction(11, 4)
GAP6_LOWER = Fraction(1, 432)
COS4_COEFF = Fraction(1, 12)
COS6_LOSS = Fraction(1, 360)
SINE_EXP_DENOM = 21


def product_sign(signs: tuple[int, ...]) -> int:
    sign = 1
    for value in signs:
        sign *= value
    return sign


def e8_positive_roots() -> list[tuple[Fraction, ...]]:
    roots: list[tuple[Fraction, ...]] = []
    zero = Fraction(0)
    one = Fraction(1)
    half = Fraction(1, 2)

    for i in range(8):
        for j in range(i + 1, 8):
            for sj in (-1, 1):
                root = [zero] * 8
                root[i] = one
                root[j] = Fraction(sj)
                roots.append(tuple(root))

    for signs in product((-1, 1), repeat=8):
        if signs[0] == 1 and product_sign(signs) == 1:
            roots.append(tuple(Fraction(sign) * half for sign in signs))

    if len(roots) != 120:
        raise AssertionError(f"expected 120 positive E8 roots, got {len(roots)}")
    return roots


def inner(a: tuple[Fraction, ...], b: tuple[Fraction, ...]) -> Fraction:
    return sum(x * y for x, y in zip(a, b))


def root_identity_check() -> bool:
    roots = e8_positive_roots()
    beta = roots[0]
    q_beta = sum(inner(alpha, beta) ** 2 for alpha in roots)
    p4_beta = sum(inner(alpha, beta) ** 4 for alpha in roots)
    p6_beta = sum(inner(alpha, beta) ** 6 for alpha in roots)
    return (
        q_beta == 60
        and p4_beta == 72
        and p6_beta == 120
        and QUARTIC_DENOM * p4_beta == q_beta * q_beta
        and SEXTIC_DENOM * p6_beta == q_beta * q_beta * q_beta
    )


def ceil_fraction(value: Fraction) -> int:
    if value < 0:
        raise ValueError("ceil_fraction expects a nonnegative value")
    return (value.numerator + value.denominator - 1) // value.denominator


def cap_value() -> Fraction:
    # E8 has kappa = 30, so S_delta(n) = delta^2 * 30n/(4d).
    return DELTA * DELTA * 30 * N / (4 * D)


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
    while s0 + GRID_STEP <= cap:
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
    exp_ceiling = ceil_fraction(s1 + t1)
    s_integral = (s1**SHAPE - s0**SHAPE) / SHAPE
    t_integral = (t1**SHAPE - t0**SHAPE) / SHAPE
    return (
        gap_lower(s0, s1, t0, t1) ** 2
        * Fraction(1, E_UPPER**exp_ceiling)
        * s_integral
        * t_integral
        * gamma_inv_square()
        / SHAPE
    )


def a0_lower() -> Fraction:
    return a0_prefactor(ROOT_DATA[GROUP]) / (PI_UPPER**RANK)


def sine_density_lower(s1: Fraction, t1: Fraction) -> Fraction:
    exponent = ceil_fraction(Fraction(4 * D, SINE_EXP_DENOM * N) * (s1 + t1))
    return Fraction(1, E_UPPER**exponent)


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
    if len(cells) != 17797:
        raise AssertionError(f"expected 17797 retained grid cells, got {len(cells)}")
    gaps = [gap_lower(s0, s1, t0, t1) for s0, s1, t0, t1 in cells]
    h_values = [local_h_lower(s0, s1, t0, t1) for s0, s1, t0, t1 in cells]
    direct_values = [local_direct_factor_lower(h) for h in h_values]
    rect_fracs = [rectangle_fraction_lower(*cell) for cell in cells]
    sine_exponents = [
        ceil_fraction(Fraction(4 * D, SINE_EXP_DENOM * N) * (s1 + t1))
        for _, s1, _, t1 in cells
    ]
    a0 = a0_lower()
    ratio131 = ratio_lower(cells, a0)
    best_term = max(positive_rectangle_lower(*cell, a0) / negative_bound() for cell in cells)

    checks = {
        "root_identity_constants": root_identity_check(),
        "cell_count_17797": len(cells) == 17797,
        "all_cells_in_cap": all(s1 < cap_value() and t1 < cap_value() for _, s1, _, t1 in cells),
        "positive_gaps": min(gaps) > 0,
        "positive_h131": min(h_values) > 0,
        "positive_direct131": min(direct_values) > 0,
        "sine_exponent_le_26": max(sine_exponents) <= 26,
        "ratio131_gt_1": ratio131 > 1,
    }

    print(f"group={GROUP} n={N} delta={DELTA}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START} t_end={T_END}")
    print(f"retained_cells={len(cells)} cap={cap_value()}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"gap131_min={min(gaps)}")
    print(f"h131_min={min(h_values)} direct_factor131_min={min(direct_values)}")
    print(f"sine_exponent_max={max(sine_exponents)}")
    print(f"log10_ratio131_lower={log10_fraction(ratio131):.2f}")
    print(f"log10_best_cell_ratio131={log10_fraction(best_term):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 n=131 rectangular certificate failed")


if __name__ == "__main__":
    main()
