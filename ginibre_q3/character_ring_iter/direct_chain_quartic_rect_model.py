#!/usr/bin/env python3
"""Quartic-improved rectangular model for E7/E8 direct Chain tails.

This diagnostic upgrades `direct_chain_rect_lower_model.py` by using the
exceptional quartic identities

    E7: sum_{alpha>0} (alpha.u)^4 = Q(u)^2/27,
    E8: sum_{alpha>0} (alpha.u)^4 = Q(u)^2/50.

The E7 result explains the fixed certificate in `DCT_RECT_TRIG_E7.md`.
The E8 result shows that the same rectangle-union method is far too small.
"""

from __future__ import annotations

from math import log10, pi

import mpmath as mp

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, alpha_value, log10_fraction
from direct_chain_rect_lower_model import sine_loss_rate


QUARTIC_DENOM = {"E7": mp.mpf(27), "E8": mp.mpf(50)}


def logsum(values: list[mp.mpf]) -> mp.mpf:
    if not values:
        return mp.mpf("-inf")
    m = max(values)
    return m + mp.log10(sum(mp.power(10, value - m) for value in values))


def gamma_moment(shape: mp.mpf, power: int, lo: mp.mpf, hi: mp.mpf) -> mp.mpf:
    return mp.gammainc(shape + power, lo, hi) / mp.gamma(shape)


def rectangle_log_piece(
    data: dict[str, object],
    n: int,
    quartic_denom: mp.mpf,
    s0: mp.mpf,
    s1: mp.mpf,
    t0: mp.mpf,
    t1: mp.mpf,
) -> mp.mpf:
    d = int(data["d"])
    shape = mp.mpf(d) / 2
    rho = mp.mpf(d) / (6 * quartic_denom * n)
    if s0 - t1 - rho * s1 * s1 <= 0:
        return mp.mpf("-inf")

    gamma4 = mp.mpf(1) / 12 - mp.pi * mp.pi / 360
    h_min = 1 - (s1 + t1) / n
    h_min += gamma4 * (2 * d / quartic_denom) * (s0 * s0 + t0 * t0) / (n * n)
    if h_min <= 0:
        return mp.mpf("-inf")

    direct_factor = (2 * d) ** 2 * h_min * h_min - 4
    if direct_factor <= 0:
        return mp.mpf("-inf")

    s_m = [gamma_moment(shape, power, s0, s1) for power in range(5)]
    t_m = [gamma_moment(shape, power, t0, t1) for power in range(3)]
    integral = (
        (s_m[2] - 2 * rho * s_m[3] + rho * rho * s_m[4]) * t_m[0]
        - 2 * (s_m[1] - rho * s_m[2]) * t_m[1]
        + s_m[0] * t_m[2]
    )
    fraction = integral / shape
    if fraction <= 0:
        return mp.mpf("-inf")

    sine_loss = mp.e ** (-(4 * sine_loss_rate() * d / n) * (s1 + t1))
    return (
        mp.log10(fraction)
        + n * mp.log10(h_min)
        + mp.log10(direct_factor)
        + mp.log10(sine_loss)
    )


def group_margin(group: str, grid: int) -> tuple[int, int, mp.mpf, list[tuple[mp.mpf, int, int]]]:
    data = ROOT_DATA[group]
    n = int(data["n_tail_start"]) - 2
    d = int(data["d"])
    c = int(data["C"])
    alpha = alpha_value(data)
    r = int(data["r"])
    cap = mp.pi * mp.pi * int(data["kappa"]) * n / (4 * d)
    rho = mp.mpf(d) / (6 * QUARTIC_DENOM[group] * n)
    common = (
        log10_fraction(a0_prefactor(data))
        - r * mp.log10(mp.pi)
        + n * mp.log10(2 * d)
        - alpha * mp.log10(n)
        - mp.log10(2 * ((2 * c) ** 2 + 4))
        - n * mp.log10(2 * c)
    )

    pieces: list[tuple[mp.mpf, int, int]] = []
    for i in range(grid):
        s0 = cap * i / grid
        s1 = cap * (i + 1) / grid
        if s0 <= 0:
            continue
        max_t = max(mp.mpf(0), s0 - rho * s1 * s1)
        jmax = int(mp.floor(max_t / cap * grid))
        for j in range(jmax):
            t0 = cap * j / grid
            t1 = cap * (j + 1) / grid
            piece = rectangle_log_piece(data, n, QUARTIC_DENOM[group], s0, s1, t0, t1)
            if piece != mp.mpf("-inf"):
                pieces.append((piece, i, j))

    margin = common + logsum([piece[0] for piece in pieces])
    pieces.sort(reverse=True, key=lambda row: row[0])
    return n, len(pieces), margin, pieces[:3]


def main() -> None:
    mp.mp.dps = 80
    grid = 80
    print(f"{'G':<4}{'nD':>5}{'grid':>7}{'cells':>8}{'margin':>12}{'top cells':>24}")
    for group in ("E7", "E8"):
        n, cells, margin, top = group_margin(group, grid)
        top_cells = ",".join(f"({i},{j})" for _, i, j in top)
        print(
            f"{group:<4}{n:>5}{grid:>7}{cells:>8}{float(margin):>12.2f}{top_cells:>24}",
            flush=True,
        )


if __name__ == "__main__":
    main()
