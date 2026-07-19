#!/usr/bin/env python3
"""Truncated Macdonald-Mehta mass diagnostic for DCT-EffectiveStart.

The full A0 coefficient uses the whole Gaussian/Macdonald-Mehta integral.
At finite n, a local saddle proof can only use a root-angle cap such as

    |alpha.v| <= delta   for every positive root alpha.

For simply-laced exceptionals with Q(v)=kappa |v|^2, the rescaled radial
variable s = Q(u)/(2d) must satisfy

    s <= delta^2 kappa n / (4d)

to guarantee the cap after v = u/sqrt(n).  This script computes the fraction
of the Gamma-Wronskian radial integral inside that cap, both for Q_3 and for
the direct Chain difference D_G(n)=Q_3(n+2)-4Q_3(n).  It is diagnostic: it
does not include the remaining finite-n Taylor factors needed for a proof.
"""

from __future__ import annotations

from math import log10, pi

import mpmath as mp

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, alpha_value, log10_fraction
from verify_chain import Q3, parse_moments


LOGS = {
    "E6": "ginibre_q3/character_ring_iter/logs/e6_80.log",
    "E7": "ginibre_q3/character_ring_iter/logs/e7_70.log",
    "E8": "ginibre_q3/character_ring_iter/logs/e8_70.log",
}

DELTAS = (1.0, 1.5, 2.0, pi)


def gamma_p(a: mp.mpf, x: mp.mpf) -> mp.mpf:
    """Regularized lower incomplete gamma P(a,x)."""

    if x <= 0:
        return mp.mpf("0")
    return mp.gammainc(a, 0, x, regularized=True)


def wronskian_fraction(a: mp.mpf, s_cap: mp.mpf) -> mp.mpf:
    """Fraction of the radial Gamma-Wronskian integral with both radii capped."""

    p0 = gamma_p(a, s_cap)
    p1 = gamma_p(a + 1, s_cap)
    p2 = gamma_p(a + 2, s_cap)
    return (a + 1) * p0 * p2 - a * p1 * p1


def q3_coeff_log(group: str, n: int, q3: int) -> float:
    data = ROOT_DATA[group]
    return log10(q3) - n * log10(2 * int(data["d"])) + alpha_value(data) * log10(n)


def d_coeff_log(group: str, n: int, diff: int) -> float:
    data = ROOT_DATA[group]
    return log10(diff) - n * log10(2 * int(data["d"])) + alpha_value(data) * log10(n)


def s_cap_for(data: dict[str, object], n: int, delta: float) -> mp.mpf:
    d = int(data["d"])
    kappa = int(data["kappa"])
    return mp.mpf(delta) * mp.mpf(delta) * kappa * n / (4 * d)


def truncated_d_fraction(data: dict[str, object], n: int, delta: float) -> mp.mpf:
    d = int(data["d"])
    alpha = alpha_value(data)
    a = mp.mpf(d) / 2
    f_n = wronskian_fraction(a, s_cap_for(data, n, delta))
    f_np2 = wronskian_fraction(a, s_cap_for(data, n + 2, delta))
    ratio = mp.mpf(n) / mp.mpf(n + 2)
    return 4 * d * d * (ratio**alpha) * f_np2 - 4 * f_n


def main() -> None:
    mp.mp.dps = 100
    print(
        f"{'G':<4}{'nQ':>5}{'nD':>5}{'delta':>8}{'S_Q':>10}"
        f"{'log10 F_Q':>12}{'Qmodel-Q':>11}{'log10 F_D':>12}{'Dmodel-D':>11}"
    )
    for group, log_path in LOGS.items():
        data = ROOT_DATA[group]
        n_q = int(data["n_tail_start"])
        n_d = n_q - 2
        d = int(data["d"])
        a = mp.mpf(d) / 2
        moms = parse_moments(log_path)
        qcoef = q3_coeff_log(group, n_q, Q3(n_q, moms))
        qd = Q3(n_d, moms)
        qd2 = Q3(n_d + 2, moms)
        diff = qd2 - 4 * qd
        if diff <= 0:
            raise AssertionError(f"{group} Chain diff is not positive at n={n_d}")
        dcoef = d_coeff_log(group, n_d, diff)
        a0_log = log10_fraction(a0_prefactor(data)) - int(data["r"]) * log10(pi)
        for delta in DELTAS:
            s_cap = s_cap_for(data, n_q, delta)
            frac = wronskian_fraction(a, s_cap)
            frac_log = float(mp.log10(frac)) if frac > 0 else float("-inf")
            q_model_log = a0_log + frac_log
            d_frac = truncated_d_fraction(data, n_d, delta)
            d_frac_log = float(mp.log10(d_frac)) if d_frac > 0 else float("-inf")
            d_model_log = a0_log + d_frac_log
            print(
                f"{group:<4}{n_q:>5}{n_d:>5}{delta:>8.3f}{float(s_cap):>10.2f}"
                f"{frac_log:>12.2f}{q_model_log - qcoef:>11.2f}"
                f"{d_frac_log:>12.2f}{d_model_log - dcoef:>11.2f}"
            )


if __name__ == "__main__":
    main()
