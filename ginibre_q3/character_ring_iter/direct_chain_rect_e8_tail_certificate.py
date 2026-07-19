#!/usr/bin/env python3
"""Exact-rational certificate for the E8 rectangular direct-Chain tail.

This checks the arithmetic part of the DCT-RectTrig(E8-tail) leaf.  The
local input is the E8 quartic root identity

    sum_{alpha>0} (alpha.u)^4 = Q(u)^2 / 50.

The positive region is a finite union of 1/40 by 1/40 rectangles:

    s0 in {37, 1481/40, ..., 793/20},
    t0 in {26, 1041/40, ..., 1239/40},

with s1=s0+1/40, t1=t0+1/40, retaining exactly the cells whose rational
gap lower bound is positive.  At n=133 this produces 13558 cells.  All
transcendental inputs are replaced by rational one-sided bounds:

    333/106 < pi < 355/113,     e < 11/4,
    1/12 - pi^2/360 > 1/18,
    sin(z/2)/(z/2) >= exp(-z^2/21) on |z| <= pi.

The script compares the resulting rectangle lower bound with the
negative-region bound 520*16^n and checks monotone odd-step propagation.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import factorial, log10

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, log10_fraction


GROUP = "E8"
N0 = 133
D = 248
C_NEG = 8
RANK = 8
ALPHA = 250
SHAPE = 124
QUARTIC_DENOM = 50

GRID_STEP = Fraction(1, 40)
S_START = Fraction(37)
T_START = Fraction(26)
T_END = Fraction(31)

PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
E_UPPER = Fraction(11, 4)
GAMMA4_LOWER = Fraction(1, 18)
SINE_EXP_DENOM = 21


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


def product_sign(signs: tuple[int, ...]) -> int:
    sign = 1
    for value in signs:
        sign *= value
    return sign


def inner(a: tuple[Fraction, ...], b: tuple[Fraction, ...]) -> Fraction:
    return sum(x * y for x, y in zip(a, b))


def quartic_identity_check() -> bool:
    roots = e8_positive_roots()
    beta = roots[0]
    q_beta = sum(inner(alpha, beta) ** 2 for alpha in roots)
    p_beta = sum(inner(alpha, beta) ** 4 for alpha in roots)
    return q_beta == 60 and p_beta == 72 and QUARTIC_DENOM * p_beta == q_beta * q_beta


def cap_lower_at_n0() -> Fraction:
    # E8 has kappa = 30, so S_pi(n) = pi^2 * 30n/(4d).
    return PI_LOWER * PI_LOWER * 30 * N0 / (4 * D)


def ceil_fraction(value: Fraction) -> int:
    if value < 0:
        raise ValueError("ceil_fraction expects a nonnegative value")
    return (value.numerator + value.denominator - 1) // value.denominator


def candidate_cells() -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    cells: list[tuple[Fraction, Fraction, Fraction, Fraction]] = []
    cap = cap_lower_at_n0()
    s0 = S_START
    while s0 + GRID_STEP <= cap:
        t0 = T_START
        while t0 + GRID_STEP <= T_END:
            s1 = s0 + GRID_STEP
            t1 = t0 + GRID_STEP
            if gap_lower(s0, s1, t0, t1, N0) > 0:
                h = local_h_lower(s0, s1, t0, t1, N0)
                if h > 0 and local_direct_factor_lower(h) > 0:
                    cells.append((s0, s1, t0, t1))
            t0 += GRID_STEP
        s0 += GRID_STEP
    return cells


def gap_lower(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, n: int) -> Fraction:
    rho = Fraction(D, 6 * QUARTIC_DENOM * n)
    return s0 - t1 - rho * s1 * s1


def local_h_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, n: int
) -> Fraction:
    quartic_gain = GAMMA4_LOWER * Fraction(2 * D, QUARTIC_DENOM)
    return Fraction(1) - (s1 + t1) / n + quartic_gain * (s0 * s0 + t0 * t0) / (
        n * n
    )


def h_monotone_from_n0(s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction) -> bool:
    quartic_gain = GAMMA4_LOWER * Fraction(2 * D, QUARTIC_DENOM)
    a_coef = s1 + t1
    b_coef = quartic_gain * (s0 * s0 + t0 * t0)
    return N0 * a_coef > 2 * b_coef


def local_direct_factor_lower(h: Fraction) -> Fraction:
    return Fraction((2 * D) ** 2) * h * h - 4


def gamma_inv_square() -> Fraction:
    # Gamma(124) = 123!.
    return Fraction(1, factorial(SHAPE - 1) ** 2)


def rectangle_fraction_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, n: int
) -> Fraction:
    exp_ceiling = ceil_fraction(s1 + t1)
    s_integral = (s1**SHAPE - s0**SHAPE) / SHAPE
    t_integral = (t1**SHAPE - t0**SHAPE) / SHAPE
    return (
        gap_lower(s0, s1, t0, t1, n) ** 2
        * Fraction(1, E_UPPER**exp_ceiling)
        * s_integral
        * t_integral
        * gamma_inv_square()
        / SHAPE
    )


def a0_lower() -> Fraction:
    return a0_prefactor(ROOT_DATA[GROUP]) / (PI_UPPER**RANK)


def sine_density_lower(s1: Fraction, t1: Fraction, n: int) -> Fraction:
    exponent = ceil_fraction(Fraction(4 * D, SINE_EXP_DENOM * n) * (s1 + t1))
    return Fraction(1, E_UPPER**exponent)


def negative_bound(n: int) -> int:
    return 2 * ((2 * C_NEG) ** 2 + 4) * (2 * C_NEG) ** n


def positive_rectangle_lower(
    s0: Fraction, s1: Fraction, t0: Fraction, t1: Fraction, n: int
) -> Fraction:
    h = local_h_lower(s0, s1, t0, t1, n)
    direct = local_direct_factor_lower(h)
    if h <= 0 or direct <= 0:
        return Fraction(0)
    return (
        a0_lower()
        * Fraction((2 * D) ** n, n**ALPHA)
        * rectangle_fraction_lower(s0, s1, t0, t1, n)
        * h**n
        * direct
        * sine_density_lower(s1, t1, n)
    )


def ratio_lower(n: int) -> Fraction:
    total = sum(
        positive_rectangle_lower(s0, s1, t0, t1, n)
        for s0, s1, t0, t1 in candidate_cells()
    )
    return total / negative_bound(n)


def odd_step_ratio_lower() -> Fraction:
    h0 = min(local_h_lower(s0, s1, t0, t1, N0) for s0, s1, t0, t1 in candidate_cells())
    base = Fraction(D, C_NEG) * h0
    return base * base * Fraction(N0, N0 + 2) ** ALPHA


def main() -> None:
    cells = candidate_cells()
    if len(cells) != 13558:
        raise AssertionError(f"expected 13558 retained grid cells, got {len(cells)}")
    h_values = [local_h_lower(s0, s1, t0, t1, N0) for s0, s1, t0, t1 in cells]
    direct_values = [local_direct_factor_lower(h) for h in h_values]
    rect_fracs = [
        rectangle_fraction_lower(s0, s1, t0, t1, N0) for s0, s1, t0, t1 in cells
    ]
    ratio0 = ratio_lower(N0)
    step = odd_step_ratio_lower()
    cap_lower = cap_lower_at_n0()
    best_term = max(
        positive_rectangle_lower(s0, s1, t0, t1, N0) / negative_bound(N0)
        for s0, s1, t0, t1 in cells
    )

    checks = {
        "quartic_identity_constant": quartic_identity_check(),
        "cell_count_13558": len(cells) == 13558,
        "all_cells_in_cap": all(s1 < cap_lower and t1 < cap_lower for _, s1, _, t1 in cells),
        "positive_gaps": all(gap_lower(s0, s1, t0, t1, N0) > 0 for s0, s1, t0, t1 in cells),
        "positive_h133": min(h_values) > 0,
        "positive_direct133": min(direct_values) > 0,
        "h_monotone_from_n133": all(h_monotone_from_n0(s0, s1, t0, t1) for s0, s1, t0, t1 in cells),
        "ratio133_gt_1": ratio0 > 1,
        "odd_step_ratio_gt_1": step > 1,
    }

    print(f"group={GROUP} n0={N0}")
    print(f"grid_step={GRID_STEP} s_start={S_START} t_start={T_START} t_end={T_END}")
    print(f"retained_cells={len(cells)}")
    print(f"cap_lower_from_pi_gt_333_106={cap_lower}")
    print(f"log10_union_fraction_lower={log10_fraction(sum(rect_fracs)):.2f}")
    print(f"h133_min={min(h_values)} direct_factor133_min={min(direct_values)}")
    print(f"log10_ratio133_lower={log10_fraction(ratio0):.2f}")
    print(f"log10_best_cell_ratio133={log10_fraction(best_term):.2f}")
    print(f"log10_odd_step_ratio_lower={log10_fraction(step):.2f}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 rectangular tail certificate failed")


if __name__ == "__main__":
    main()
