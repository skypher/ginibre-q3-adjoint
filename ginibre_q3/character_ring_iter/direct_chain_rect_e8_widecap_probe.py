#!/usr/bin/env python3
"""Wide-cap E8 next-boundary diagnostic for DCT-RectTrig.

This is a proof-search diagnostic companion to the exact n=131 certificate.
It tests the same wide-cap region after the maintained n>=133 certificate by
replacing the strict pi radial cap with a slightly wider root-angle cap
delta=19/6 and by testing the finite-n constants suggested by the exact
arithmetic:

    2(1-cos z) >= z^2 - z^4/12 + z^6/432,
    2 cos z >= 2 - z^2 + z^4/12 - z^6/360,
    sin(z/2)/(z/2) >= exp(-z^2/21),
    exp(-x) >= (1 - x/M)^M.

The last line is evaluated in log-space for M in {2048,4096}.  The
exact-rational replay using e < 11/4 is
direct_chain_rect_e8_n131_certificate.py.
"""

from __future__ import annotations

from fractions import Fraction
from math import exp, lgamma, log, log10

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
GAP6 = Fraction(1, 432)
COS4_COEFF = Fraction(1, 12)
COS6_LOSS = Fraction(1, 360)
SINE_DENOM = 21


def logdiffpow(hi: float, lo: float, power: int) -> float:
    log_hi = power * log(hi)
    log_lo = power * log(lo)
    if log_lo >= log_hi:
        return float("-inf")
    return (log_hi + log(1 - exp(log_lo - log_hi))) / log(10)


def logsum10(values: list[float]) -> float:
    if not values:
        return float("-inf")
    top = max(values)
    return top + log10(sum(10 ** (value - top) for value in values))


def exp_lower_log10(x: float, m: int) -> float:
    if x >= m:
        return float("-inf")
    return m * log10(1 - x / m)


def cap_value() -> Fraction:
    return DELTA * DELTA * 30 * N / (4 * D)


def cell_log_margin(
    s0_q: Fraction,
    s1_q: Fraction,
    t0_q: Fraction,
    t1_q: Fraction,
    exp_m: int,
) -> float:
    s0 = float(s0_q)
    s1 = float(s1_q)
    t0 = float(t0_q)
    t1 = float(t1_q)

    rho = D / (6 * QUARTIC_DENOM * N)
    gap = s0 - t1 - rho * s1 * s1
    gap += float(GAP6) * ((2 * D) ** 2) * s0**3 / (SEXTIC_DENOM * N * N)
    if gap <= 0:
        return float("-inf")

    h = 1 - (s1 + t1) / N
    h += float(COS4_COEFF) * (2 * D / QUARTIC_DENOM) * (s0 * s0 + t0 * t0) / (N * N)
    h -= float(COS6_LOSS) * ((2 * D) ** 2 / SEXTIC_DENOM) * (s1**3 + t1**3) / (N**3)
    if h <= 0:
        return float("-inf")
    direct = (2 * D) ** 2 * h * h - 4
    if direct <= 0:
        return float("-inf")

    rect_log = 2 * log10(gap)
    rect_log += exp_lower_log10(s1 + t1, exp_m)
    rect_log += logdiffpow(s1, s0, SHAPE) - log10(SHAPE)
    rect_log += logdiffpow(t1, t0, SHAPE) - log10(SHAPE)
    rect_log -= 2 * lgamma(SHAPE) / log(10)
    rect_log -= log10(SHAPE)

    sine_arg = (4 * D / (SINE_DENOM * N)) * (s1 + t1)
    a0_log = log10_fraction(a0_prefactor(ROOT_DATA[GROUP])) - RANK * log10(float(PI_UPPER))
    negative_log = log10(2 * ((2 * C_NEG) ** 2 + 4)) + N * log10(2 * C_NEG)
    return (
        a0_log
        + N * log10(2 * D)
        - ALPHA * log10(N)
        + rect_log
        + N * log10(h)
        + log10(direct)
        + exp_lower_log10(sine_arg, exp_m)
        - negative_log
    )


def scan(exp_m: int) -> tuple[int, float, float]:
    cap = cap_value()
    values: list[float] = []
    s0 = S_START
    while s0 + GRID_STEP <= cap:
        t0 = T_START
        while t0 + GRID_STEP <= T_END:
            value = cell_log_margin(s0, s0 + GRID_STEP, t0, t0 + GRID_STEP, exp_m)
            if value != float("-inf"):
                values.append(value)
            t0 += GRID_STEP
        s0 += GRID_STEP
    return len(values), max(values), logsum10(values)


def main() -> None:
    print(f"group={GROUP} n={N} delta={DELTA} grid_step={GRID_STEP}", flush=True)
    print(f"s_start={S_START} t_start={T_START} t_end={T_END}", flush=True)
    print(f"cap={cap_value()} log10_A0_lower={log10_fraction(a0_prefactor(ROOT_DATA[GROUP]) / (PI_UPPER**RANK)):.2f}", flush=True)
    for exp_m in (2048, 4096):
        cells, best, total = scan(exp_m)
        print(
            f"exp_M={exp_m} cells={cells} "
            f"log10_best_cell={best:.2f} log10_total={total:.2f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
