#!/usr/bin/env python3
"""Exact E8 finite-bridge replay from the arXiv m_0..m_100 moment table.

The local character-ring log supplies m_0..m_70.  The arXiv 2412.21189
ancillary table supplies m_71..m_100.  This checker verifies the source
overlap through the n=81 certificate and then evaluates Q_3(n) and the direct
Chain differences

    D_E8(n) = Q_3^{E8}(n+2) - 4 Q_3^{E8}(n)

for the remaining finite E8 bridge window.
"""

from __future__ import annotations

from math import comb

import direct_chain_rect_e8_n81_certificate as moments
from direct_chain_a0_dominance import log10_fraction


MOMENTS_0_100 = moments.MOMENTS_0_100
BRIDGE_ODD_N = tuple(range(69, 78, 2))
PREFIX_CHAIN_ODD_N = tuple(range(1, 66, 2))
CHAIN_ODD_N = tuple(range(67, 78, 2))
EXTENDED_CHAIN_ODD_N = tuple(range(67, 96, 2))


def q3(n: int) -> int:
    return 2 * sum(
        comb(n, k)
        * (
            MOMENTS_0_100[k + 2] * MOMENTS_0_100[n - k]
            - MOMENTS_0_100[k + 1] * MOMENTS_0_100[n - k + 1]
        )
        for k in range(n + 1)
    )


def chain_diff(n: int) -> int:
    return q3(n + 2) - 4 * q3(n)


def main() -> None:
    prefix_chain_values = {n: chain_diff(n) for n in PREFIX_CHAIN_ODD_N}
    bridge_values = {n: q3(n) for n in BRIDGE_ODD_N}
    chain_values = {n: chain_diff(n) for n in CHAIN_ODD_N}
    extended_chain_values = {n: chain_diff(n) for n in EXTENDED_CHAIN_ODD_N}

    checks = {
        "moment_log_0_100_arxiv_extension": moments.moment_source_check(),
        "prefix_chain_diff_1_65_positive": all(
            value > 0 for value in prefix_chain_values.values()
        ),
        "bridge_q3_69_77_positive": all(value > 0 for value in bridge_values.values()),
        "chain_diff_67_77_positive": all(value > 0 for value in chain_values.values()),
        "extended_chain_diff_67_95_positive": all(
            value > 0 for value in extended_chain_values.values()
        ),
    }

    print("group=E8 source=m_0..m_100")
    print(f"prefix_chain_odd_n={PREFIX_CHAIN_ODD_N}")
    print(f"bridge_odd_n={BRIDGE_ODD_N}")
    print(f"chain_odd_n={CHAIN_ODD_N}")
    print(f"extended_chain_odd_n={EXTENDED_CHAIN_ODD_N}")
    print(f"min_prefix_chain_diff={min(prefix_chain_values.values())}")
    print(
        "log10_min_prefix_chain_diff="
        f"{min(log10_fraction(value) for value in prefix_chain_values.values()):.2f}"
    )
    print(f"log10_min_bridge_q3={min(log10_fraction(value) for value in bridge_values.values()):.2f}")
    print(f"log10_min_chain_diff={min(log10_fraction(value) for value in chain_values.values()):.2f}")
    print(
        "log10_min_extended_chain_diff="
        f"{min(log10_fraction(value) for value in extended_chain_values.values()):.2f}"
    )
    for n, value in bridge_values.items():
        print(f"Q3_{n}={value}")
    for n, value in chain_values.items():
        print(f"D_E8_{n}={value}")
    for name, ok in checks.items():
        print(f"{name}: {'OK' if ok else 'FAIL'}")
    if not all(checks.values()):
        raise AssertionError("E8 arXiv m100 finite bridge replay failed")


if __name__ == "__main__":
    main()
