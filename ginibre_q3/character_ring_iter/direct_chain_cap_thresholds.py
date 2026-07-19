#!/usr/bin/env python3
"""Solve the cap size needed by the finite-n identity-saddle model.

The truncated-mass diagnostic asks what fraction of the Macdonald-Mehta
identity saddle is available under a root-angle cap |alpha.v| <= delta.
This script inverts that comparison at the current exact prefixes.  It
finds the smallest delta for which the truncated leading model reaches the
exact logged coefficient of Q_3(n), and separately of the direct Chain
difference

    D_G(n) = Q_3(n+2) - 4 Q_3(n).

The output is diagnostic.  It does not include a finite-n trigonometric
lower bound; it identifies how wide such a lower-bound region must be.
"""

from __future__ import annotations

from math import log10, pi

import mpmath as mp

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, alpha_value, log10_fraction
from direct_chain_truncated_mass import LOGS, d_coeff_log, truncated_d_fraction
from direct_chain_truncated_mass import q3_coeff_log, s_cap_for, wronskian_fraction
from verify_chain import Q3, parse_moments


def log10_a0(data: dict[str, object]) -> mp.mpf:
    r = int(data["r"])
    return mp.mpf(log10_fraction(a0_prefactor(data))) - r * mp.log10(mp.pi)


def find_delta_for_log_target(
    log_target: mp.mpf,
    value_log_at_delta,
    hi_start: mp.mpf = mp.mpf("1"),
) -> mp.mpf:
    lo = mp.mpf("0")
    hi = hi_start
    while value_log_at_delta(hi) < log_target:
        hi *= 2
        if hi > 100:
            raise RuntimeError("failed to bracket delta threshold")

    for _ in range(120):
        mid = (lo + hi) / 2
        if value_log_at_delta(mid) < log_target:
            lo = mid
        else:
            hi = mid
    return hi


def q_fraction_log(data: dict[str, object], n: int, delta: mp.mpf) -> mp.mpf:
    a = mp.mpf(int(data["d"])) / 2
    frac = wronskian_fraction(a, s_cap_for(data, n, float(delta)))
    if frac <= 0:
        return mp.mpf("-inf")
    return mp.log10(frac)


def d_fraction_log(data: dict[str, object], n: int, delta: mp.mpf) -> mp.mpf:
    frac = truncated_d_fraction(data, n, float(delta))
    if frac <= 0:
        return mp.mpf("-inf")
    return mp.log10(frac)


def main() -> None:
    mp.mp.dps = 120
    print(
        f"{'G':<4}{'nQ':>5}{'nD':>5}{'log10 FQ*':>12}{'deltaQ':>10}"
        f"{'deltaQ/pi':>11}{'S_Q':>10}{'log10 FD*':>12}{'deltaD':>10}"
        f"{'deltaD/pi':>11}{'S_D':>10}"
    )
    for group, log_path in LOGS.items():
        data = ROOT_DATA[group]
        n_q = int(data["n_tail_start"])
        n_d = n_q - 2
        moms = parse_moments(log_path)

        q_log_target = mp.mpf(q3_coeff_log(group, n_q, Q3(n_q, moms))) - log10_a0(data)
        q_delta = find_delta_for_log_target(
            q_log_target, lambda delta: q_fraction_log(data, n_q, delta)
        )

        diff = Q3(n_d + 2, moms) - 4 * Q3(n_d, moms)
        if diff <= 0:
            raise AssertionError(f"{group} Chain diff is not positive at n={n_d}")
        d_log_target = mp.mpf(d_coeff_log(group, n_d, diff)) - log10_a0(data)
        d_delta = find_delta_for_log_target(
            d_log_target, lambda delta: d_fraction_log(data, n_d, delta)
        )

        s_q = s_cap_for(data, n_q, float(q_delta))
        s_d = s_cap_for(data, n_d, float(d_delta))
        print(
            f"{group:<4}{n_q:>5}{n_d:>5}{float(q_log_target):>12.2f}"
            f"{float(q_delta):>10.4f}{float(q_delta / pi):>11.3f}"
            f"{float(s_q):>10.2f}{float(d_log_target):>12.2f}"
            f"{float(d_delta):>10.4f}{float(d_delta / pi):>11.3f}"
            f"{float(s_d):>10.2f}"
        )


if __name__ == "__main__":
    main()
