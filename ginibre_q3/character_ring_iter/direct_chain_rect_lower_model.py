#!/usr/bin/env python3
"""Rectangular radial lower-bound model for DCT-WideRegionTrig.

This diagnostic tests a proof-shaped local estimate for

    D_G(n) = Q_3(n+2) - 4 Q_3(n)

at the current E_6/E_7/E_8 prefixes.  On a root cap |alpha.v| <= pi it uses
only elementary inequalities:

    2(1 - cos t) >= (4/pi^2) t^2,
    2(1 - cos t) <= t^2,
    2 cos t >= 2 - t^2 + (1/12 - pi^2/360)t^4,

and an exponential lower model for the Weyl sine factor.  The radial region
is a rectangle s in [s0,s1], t in [0,t1] with t1 < (4/pi^2)s0, plus its
swap.  Gamma moments on the rectangle are evaluated directly.

This is not the final theorem: it is a size test for a conservative
rectangular subregion of the desired wide-region trigonometric estimate.
"""

from __future__ import annotations

from math import log10, pi

import mpmath as mp

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, alpha_value, log10_fraction
from direct_chain_truncated_mass import LOGS, d_coeff_log, s_cap_for
from verify_chain import Q3, parse_moments


def gamma_moment(shape: mp.mpf, power: int, lo: mp.mpf, hi: mp.mpf) -> mp.mpf:
    return mp.gammainc(shape + power, lo, hi) / mp.gamma(shape)


def log10_a0(data: dict[str, object]) -> mp.mpf:
    return mp.mpf(log10_fraction(a0_prefactor(data))) - int(data["r"]) * mp.log10(mp.pi)


def sine_loss_rate() -> mp.mpf:
    """Endpoint exponential rate for sin(t/2)/(t/2) on |t| <= pi."""

    return -mp.log(mp.mpf(2) / mp.pi) / (mp.pi * mp.pi)


def rectangle_fraction(
    shape: mp.mpf,
    mu: mp.mpf,
    s0: mp.mpf,
    s1: mp.mpf,
    t1: mp.mpf,
) -> mp.mpf:
    s_m0 = gamma_moment(shape, 0, s0, s1)
    s_m1 = gamma_moment(shape, 1, s0, s1)
    s_m2 = gamma_moment(shape, 2, s0, s1)
    t_m0 = gamma_moment(shape, 0, 0, t1)
    t_m1 = gamma_moment(shape, 1, 0, t1)
    t_m2 = gamma_moment(shape, 2, 0, t1)
    wedge = mu * mu * s_m2 * t_m0 - 2 * mu * s_m1 * t_m1 + s_m0 * t_m2
    # Divide by 2a, the full leading Wronskian radial integral, and multiply
    # by 2 for the swapped rectangle.
    return wedge / shape


def rectangle_log_margin(
    data: dict[str, object],
    n: int,
    s0: mp.mpf,
    s1: mp.mpf,
    t1: mp.mpf,
) -> mp.mpf:
    d = int(data["d"])
    r_plus = int(data["R_plus"])
    shape = mp.mpf(d) / 2
    alpha = alpha_value(data)
    c = int(data["C"])
    neg_factor = 2 * ((2 * c) ** 2 + 4)
    mu = 4 / (mp.pi * mp.pi)
    gamma4 = mp.mpf(1) / 12 - mp.pi * mp.pi / 360
    quartic_gain = 2 * gamma4 * d / r_plus
    h_min = 1 - (s1 + t1) / n + quartic_gain * s0 * s0 / (n * n)
    if h_min <= 0:
        return mp.mpf("-inf")
    direct_factor = (2 * d) ** 2 * h_min * h_min - 4
    if direct_factor <= 0:
        return mp.mpf("-inf")
    frac = rectangle_fraction(shape, mu, s0, s1, t1)
    if frac <= 0:
        return mp.mpf("-inf")

    sine_loss = mp.e ** (-(4 * sine_loss_rate() * d / n) * (s1 + t1))
    lower_log = (
        log10_a0(data)
        + n * mp.log10(2 * d)
        - alpha * mp.log10(n)
        + mp.log10(frac)
        + n * mp.log10(h_min)
        + mp.log10(direct_factor)
        + mp.log10(sine_loss)
    )
    negative_log = mp.log10(neg_factor) + n * mp.log10(2 * c)
    return lower_log - negative_log


def search_group(group: str) -> tuple[int, mp.mpf, tuple[mp.mpf, mp.mpf, mp.mpf]]:
    data = ROOT_DATA[group]
    n = int(data["n_tail_start"]) - 2
    s_cap = s_cap_for(data, n, pi)
    mu = 4 / (mp.pi * mp.pi)
    best = (mp.mpf("-inf"), (mp.mpf("0"), mp.mpf("0"), mp.mpf("0")))
    for i in range(12, 75):
        s0 = s_cap * i / 80
        for width_i in range(2, 17):
            s1 = s0 + s_cap * width_i / 160
            if s1 >= s_cap:
                continue
            for tau_i in range(2, 15):
                t1 = mu * s0 * tau_i / 16
                if t1 >= s_cap:
                    continue
                margin = rectangle_log_margin(data, n, s0, s1, t1)
                if margin > best[0]:
                    best = (margin, (s0, s1, t1))
    return n, best[0], best[1]


def main() -> None:
    mp.mp.dps = 100
    print(
        f"{'G':<4}{'nD':>5}{'cap S':>10}{'best margin':>14}"
        f"{'s0':>10}{'s1':>10}{'t1':>10}{'exact D margin':>16}",
        flush=True,
    )
    for group, log_path in LOGS.items():
        data = ROOT_DATA[group]
        n, margin, box = search_group(group)
        moms = parse_moments(log_path)
        diff = Q3(n + 2, moms) - 4 * Q3(n, moms)
        exact_log = d_coeff_log(group, n, diff)
        alpha = alpha_value(data)
        d = int(data["d"])
        c = int(data["C"])
        neg_factor = 2 * ((2 * c) ** 2 + 4)
        exact_margin = (
            exact_log
            + n * log10(2 * d)
            - alpha * log10(n)
            - log10(neg_factor)
            - n * log10(2 * c)
        )
        print(
            f"{group:<4}{n:>5}{float(s_cap_for(data, n, pi)):>10.2f}"
            f"{float(margin):>14.2f}{float(box[0]):>10.2f}"
            f"{float(box[1]):>10.2f}{float(box[2]):>10.2f}"
            f"{exact_margin:>16.2f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
