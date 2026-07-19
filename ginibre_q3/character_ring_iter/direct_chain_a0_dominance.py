#!/usr/bin/env python3
"""Exact A0 constants for the exceptional direct Chain-tail test.

For the simply-laced exceptional groups E6, E7, and E8, this script computes
the leading positive-saddle coefficient in the form

    A0(G) = A0_prefactor(G) / pi^r.

The prefactor includes the normalized-torus lattice factor

    nu_G = |P^vee/Q^vee|^2 / covol(Q^vee)^2.

For the simply-laced adjoint forms below this is the Cartan determinant:
`3, 2, 1` for `E6, E7, E8`.

The computation uses the Macdonald-Mehta Gaussian integral at k = 1 with
the standard long-root normalization |alpha|^2 = 2; see MIT OCW 18.735,
Chapter 4, Theorem 4.1.  It then checks the following exact consequence of
a one-term positive-saddle estimate:

    |E_G(n)| <= A0(G)/2 and |E_G(n+2)| <= A0(G)/2

imply

    D_{G,+}(n) >= (2d)^n n^(-alpha) A0(G)
                 * (2 d^2 (n/(n+2))^alpha - 6).

The printed dominance margin uses pi < 4, so it is a rational lower bound
for A0(G).  This is a conditional checker: it shows that a half-leading
estimate would dominate the bad-region bound, but it does not prove that
the half-leading estimate holds at the current exact prefixes.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial, log10


ROOT_DATA = {
    "E6": {
        "r": 6,
        "R_plus": 36,
        "d": 78,
        "C": 6,
        "kappa": 12,
        "degrees": (2, 5, 6, 8, 9, 12),
        "weyl_lattice_factor": 3,
        "n_tail_start": 77,
    },
    "E7": {
        "r": 7,
        "R_plus": 63,
        "d": 133,
        "C": 7,
        "kappa": 18,
        "degrees": (2, 6, 8, 10, 12, 14, 18),
        "weyl_lattice_factor": 2,
        "n_tail_start": 67,
    },
    "E8": {
        "r": 8,
        "R_plus": 120,
        "d": 248,
        "C": 8,
        "kappa": 30,
        "degrees": (2, 8, 12, 14, 18, 20, 24, 30),
        "weyl_lattice_factor": 1,
        "n_tail_start": 67,
    },
}


def alpha_value(data: dict[str, object]) -> int:
    return 2 + 2 * int(data["R_plus"]) + int(data["r"])


def factorial_product_minus_one(degrees: tuple[int, ...]) -> int:
    out = 1
    for degree in degrees:
        out *= factorial(degree - 1)
    return out


def a0_prefactor(data: dict[str, object]) -> Fraction:
    """Return P with A0 = P / pi^r.

    Macdonald-Mehta gives

        Z_lambda = (d/kappa)^a (2 pi)^(r/2) prod_i d_i!

    for a = |R_+| + r/2 = d/2 and lambda = kappa/(2d).  Since
    |W| = prod_i d_i and J = 8 a d^2 Z_lambda^2, the Weyl-normalized
    leading coefficient is

        A0 = nu_G 8 a d^2 (d/kappa)^d
             prod_i((d_i-1)!)^2 / (2 pi)^r.
    """

    r = int(data["r"])
    d = int(data["d"])
    r_plus = int(data["R_plus"])
    kappa = int(data["kappa"])
    degrees = data["degrees"]
    lattice_factor = int(data["weyl_lattice_factor"])
    if not isinstance(degrees, tuple):
        raise TypeError("degrees must be a tuple")
    if Fraction(2 * r_plus, r) != kappa:
        raise AssertionError("kappa does not match trace of positive roots")
    a = Fraction(d, 2)
    prod = factorial_product_minus_one(degrees)
    return (
        Fraction(lattice_factor)
        * 8
        * a
        * d
        * d
        * Fraction(d, kappa) ** d
        * prod
        * prod
        / (2**r)
    )


def log10_fraction(value: Fraction) -> float:
    return log10(value.numerator) - log10(value.denominator)


def half_lead_factor(data: dict[str, object], n: int) -> Fraction:
    d = int(data["d"])
    alpha = alpha_value(data)
    ratio = Fraction(n, n + 2) ** alpha
    return 2 * d * d * ratio - 6


def dominance_log10_margin(data: dict[str, object], n: int) -> float:
    """Return log10(lower positive / negative bound) under half-leading input."""

    r = int(data["r"])
    d = int(data["d"])
    c = int(data["C"])
    alpha = alpha_value(data)
    neg_factor = 2 * ((2 * c) ** 2 + 4)
    lead = half_lead_factor(data, n)
    if lead <= 0:
        return float("-inf")

    a0_lower = a0_prefactor(data) / (4**r)
    return (
        log10_fraction(a0_lower)
        + log10_fraction(lead)
        + n * log10(Fraction(d, c))
        - alpha * log10(n)
        - log10(neg_factor)
    )


def first_odd_dominance(data: dict[str, object], limit: int = 999) -> int | None:
    for n in range(1, limit + 1, 2):
        if dominance_log10_margin(data, n) >= 0:
            return n
    return None


def main() -> None:
    print(
        f"{'G':<4}{'alpha':>8}{'log10 A0(pi<4)':>18}"
        f"{'first cond':>12}{'tail n':>10}{'cond margin':>14}"
    )
    for group in ("E6", "E7", "E8"):
        data = ROOT_DATA[group]
        r = int(data["r"])
        n_tail = int(data["n_tail_start"])
        if n_tail % 2 == 0:
            raise AssertionError("tail start must be odd")
        a0_lower_log = log10_fraction(a0_prefactor(data) / (4**r))
        first = first_odd_dominance(data)
        margin = dominance_log10_margin(data, n_tail)
        print(
            f"{group:<4}{alpha_value(data):>8}{a0_lower_log:>18.2f}"
            f"{str(first):>12}{n_tail:>10}{margin:>14.2f}"
        )


if __name__ == "__main__":
    main()
