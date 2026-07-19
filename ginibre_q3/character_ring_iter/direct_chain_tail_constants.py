#!/usr/bin/env python3
"""Exact constants for the direct Chain-tail reduction.

Let n be odd and let

    D_G(n) := Q_3^G(n+2) - 4 Q_3^G(n).

If the positive-saddle part of Q_3 has the local expansion

    Q_{3,+}(n) = (2d)^n n^{-alpha} (A0 + A1/n + O(n^{-2})),

then the positive-saddle part of D_G has expansion

    D_+(n) = (2d)^n n^{-alpha}
             (4(d^2-1) A0 + (-8 alpha d^2 A0 + 4(d^2-1) A1)/n
              + O(n^{-2})).

This script prints the exact multiplier and the same-base negative-region
constant inherited from the existing |I_-| bound.  It is a checker for the
reduction target, not a proof of the remaining Wong remainder bound.
"""

from __future__ import annotations

from fractions import Fraction


ROOT_DATA = {
    "G2": {"r": 2, "R_plus": 6, "d": 14, "C": 2, "nmax": 197},
    "F4": {"r": 4, "R_plus": 24, "d": 52, "C": 4, "nmax": 217},
    "E6": {"r": 6, "R_plus": 36, "d": 78, "C": 2, "nmax": 77},
    "E7": {"r": 7, "R_plus": 63, "d": 133, "C": 3, "nmax": 67},
    "E8": {"r": 8, "R_plus": 120, "d": 248, "C": 8, "nmax": 67},
}


def alpha_value(r: int, r_plus: int) -> int:
    return 2 + 2 * r_plus + r


def constants(group: str) -> dict[str, Fraction | int]:
    rd = ROOT_DATA[group]
    d = rd["d"]
    c = rd["C"]
    alpha = alpha_value(rd["r"], rd["R_plus"])
    leading_multiplier = 4 * (d * d - 1)
    first_geom = -8 * alpha * d * d
    inherited_negative_factor = 2 * ((2 * c) ** 2 + 4)
    cartan_ratio = Fraction(d, c)
    return {
        "alpha": alpha,
        "leading_multiplier": leading_multiplier,
        "first_geom": first_geom,
        "negative_factor": inherited_negative_factor,
        "cartan_ratio": cartan_ratio,
        "nmax": rd["nmax"],
    }


def main() -> None:
    print(
        f"{'G':<4}{'alpha':>8}{'4(d^2-1)':>14}{'-8 alpha d^2':>18}"
        f"{'neg factor':>14}{'d/C':>12}{'nmax':>8}"
    )
    for group in ["G2", "F4", "E6", "E7", "E8"]:
        c = constants(group)
        print(
            f"{group:<4}{c['alpha']:>8}{c['leading_multiplier']:>14}"
            f"{c['first_geom']:>18}{c['negative_factor']:>14}"
            f"{str(c['cartan_ratio']):>12}{c['nmax']:>8}"
        )


if __name__ == "__main__":
    main()
