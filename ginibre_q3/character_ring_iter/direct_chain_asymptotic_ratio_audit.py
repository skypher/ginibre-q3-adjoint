#!/usr/bin/env python3
"""Exact-log diagnostics for the direct Chain asymptotic scale.

This script compares exact logged values of Q_3(n) and

    D_G(n) = Q_3(n+2) - 4 Q_3(n)

against the Macdonald-Mehta leading coefficient A0(G).  It does not prove a
tail estimate.  Its purpose is to locate the effective-start obstruction:
at the current logged n, the asymptotic coefficient has not yet emerged.
"""

from __future__ import annotations

from math import log10

from direct_chain_a0_dominance import (
    ROOT_DATA,
    a0_prefactor,
    alpha_value,
    half_lead_factor,
    log10_fraction,
)
from verify_chain import Q3, parse_moments


LOGS = {
    "E6": "ginibre_q3/character_ring_iter/logs/e6_80.log",
    "E7": "ginibre_q3/character_ring_iter/logs/e7_70.log",
    "E8": "ginibre_q3/character_ring_iter/logs/e8_70.log",
}


def log10_int(value: int) -> float:
    if value <= 0:
        raise ValueError("positive integer expected")
    return log10(value)


def q3_coeff_log(group: str, n: int, q3: int) -> float:
    data = ROOT_DATA[group]
    d = int(data["d"])
    alpha = alpha_value(data)
    return log10_int(q3) - n * log10(2 * d) + alpha * log10(n)


def d_coeff_log(group: str, n: int, diff: int) -> float:
    data = ROOT_DATA[group]
    d = int(data["d"])
    alpha = alpha_value(data)
    return log10_int(diff) - n * log10(2 * d) + alpha * log10(n)


def main() -> None:
    print(
        f"{'G':<4}{'nQ':>6}{'log10 A0<=':>14}{'log10 Qcoef':>14}"
        f"{'A0-Qcoef':>12}{'nD':>6}{'log10 Dcoef':>14}{'leadD-Dcoef':>14}"
    )
    for group, log_path in LOGS.items():
        data = ROOT_DATA[group]
        n_q = int(data["n_tail_start"])
        moms = parse_moments(log_path)
        n_d = len(moms) - 5
        if n_d % 2 == 0:
            n_d -= 1
        qn = Q3(n_q, moms)
        qd = Q3(n_d, moms)
        qd2 = Q3(n_d + 2, moms)
        diff = qd2 - 4 * qd
        if diff <= 0:
            raise AssertionError(f"{group} Chain diff is not positive at n={n_d}")

        r = int(data["r"])
        a0_log_upper = log10_fraction(a0_prefactor(data))
        # A0 = prefactor / pi^r, so pi > 3 gives a simple upper estimate.
        a0_log_upper -= r * log10(3)
        qcoef_log = q3_coeff_log(group, n_q, qn)
        dcoef_log = d_coeff_log(group, n_d, diff)
        lead_d_log = a0_log_upper + log10_fraction(half_lead_factor(data, n_d))
        print(
            f"{group:<4}{n_q:>6}{a0_log_upper:>14.2f}{qcoef_log:>14.2f}"
            f"{a0_log_upper - qcoef_log:>12.2f}{n_d:>6}{dcoef_log:>14.2f}"
            f"{lead_d_log - dcoef_log:>14.2f}"
        )


if __name__ == "__main__":
    main()
