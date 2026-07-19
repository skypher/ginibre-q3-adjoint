#!/usr/bin/env python3
"""Audit the half-leading direct Chain-tail target against exact logs.

For odd n, write Q_3(n) = Q_{3,+}(n) - |I_-(n)| with

    |I_-(n)| <= 2(2C_G)^n.

The proposed prefix-level half-leading estimate would imply

    Q_{3,+}(n) >= (2d)^n n^(-alpha) A0(G)/2.

Since Q_{3,+}(n) <= Q_3(n) + 2(2C_G)^n, the estimate is impossible at a
tested n if the lower bound exceeds this exact upper bound.  This script
checks that comparison using exact rational arithmetic wherever a
character-ring log supplies the needed moments.
"""

from __future__ import annotations

from fractions import Fraction
from math import log10

from direct_chain_a0_dominance import ROOT_DATA, a0_prefactor, alpha_value, log10_fraction
from verify_chain import Q3, parse_moments


LOGS = {
    "E6": "ginibre_q3/character_ring_iter/logs/e6_80.log",
    "E7": "ginibre_q3/character_ring_iter/logs/e7_70.log",
    "E8": "ginibre_q3/character_ring_iter/logs/e8_70.log",
}


def log10_int(value: int) -> float:
    if value <= 0:
        raise ValueError("log10_int expects a positive integer")
    return log10(value)


def halflead_lower_bound(data: dict[str, object], n: int) -> Fraction:
    d = int(data["d"])
    r = int(data["r"])
    alpha = alpha_value(data)
    return a0_prefactor(data) * (2 * d) ** n / (2 * 4**r * n**alpha)


def audit_group(group: str, log_path: str) -> dict[str, object]:
    data = ROOT_DATA[group]
    c = int(data["C"])
    moms = parse_moments(log_path)
    last_available_odd = len(moms) - 4
    if last_available_odd % 2 == 0:
        last_available_odd -= 1

    rows = []
    impossible = []
    not_ruled = []
    for n in range(1, last_available_odd + 1, 2):
        q3 = Q3(n, moms)
        i_minus_bound = 2 * (2 * c) ** n
        qplus_upper = q3 + i_minus_bound
        halflead = halflead_lower_bound(data, n)
        is_impossible = halflead > qplus_upper
        excess = log10_fraction(halflead) - log10_int(qplus_upper)
        row = (n, qplus_upper, halflead, excess, is_impossible)
        rows.append(row)
        if is_impossible:
            impossible.append(row)
        else:
            not_ruled.append(row)
    return {
        "rows": rows,
        "impossible": impossible,
        "not_ruled": not_ruled,
        "last_available_odd": last_available_odd,
    }


def main() -> None:
    print(
        f"{'G':<4}{'n':>6}{'log10 Q3+Ibd':>16}"
        f"{'log10 halflead':>18}{'excess':>12}{'status':>14}"
    )
    for group, log_path in LOGS.items():
        result = audit_group(group, log_path)
        n = int(ROOT_DATA[group]["n_tail_start"])
        row_by_n = {row[0]: row for row in result["rows"]}
        _, qplus_upper, halflead, excess, is_impossible = row_by_n[n]
        status = "IMPOSSIBLE" if is_impossible else "not ruled"
        print(
            f"{group:<4}{n:>6}{log10_int(qplus_upper):>16.2f}"
            f"{log10_fraction(halflead):>18.2f}{excess:>12.2f}{status:>14}"
        )

    print("\nFull exact-log scan:")
    print(
        f"{'G':<4}{'tested odds':>14}{'last impossible':>18}"
        f"{'first not ruled':>18}{'min excess':>14}"
    )
    for group, log_path in LOGS.items():
        result = audit_group(group, log_path)
        impossible = result["impossible"]
        not_ruled = result["not_ruled"]
        last_impossible = impossible[-1][0] if impossible else None
        first_not_ruled = not_ruled[0][0] if not_ruled else None
        min_excess = min(row[3] for row in result["rows"])
        print(
            f"{group:<4}{'1..' + str(result['last_available_odd']):>14}"
            f"{str(last_impossible):>18}{str(first_not_ruled):>18}"
            f"{min_excess:>14.2f}"
        )


if __name__ == "__main__":
    main()
