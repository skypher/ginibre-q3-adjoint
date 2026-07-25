#!/usr/bin/env python3
"""Exact limit of the difference-factorization reduction.

Two recent results use the same move:

    V_1^2 = V_0 + V_2          =>  D_2 = D_1 S_1        (ordinary ring)
    T^2   = B_0 + B_1          =>  D_1 = D_T S_T        (odd orbit ring)

Both convert a minus factor into one minus and one plus factor of a smaller
label, preserving minus parity, uniformly in level and in word length.

This program verifies that in the *ordinary* representation ring the move has
exactly one instance.  Everything is exact integer bivariate polynomial
arithmetic; there is no floating point and no external dependency.

Checks:
  1. Clebsch-Gordan squares:  chi_u^2 = sum_(i=0..u) chi_2i, with u+1 summands.
  2. The known instance:      D_2 = D_1 S_1 exactly.
  3. Uniqueness:              D_a != D_u S_u for every other (a,u) in range.
"""

from __future__ import annotations

import sys

# Bivariate integer polynomials as {(i, j): coeff} for x^i y^j.
Poly = dict


def padd(p: Poly, q: Poly, s: int = 1) -> Poly:
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + s * c
        if r[m] == 0:
            del r[m]
    return r


def pmul(p: Poly, q: Poly) -> Poly:
    r: Poly = {}
    for (i1, j1), c1 in p.items():
        for (i2, j2), c2 in q.items():
            m = (i1 + i2, j1 + j2)
            r[m] = r.get(m, 0) + c1 * c2
            if r[m] == 0:
                del r[m]
    return r


ONE: Poly = {(0, 0): 1}
X: Poly = {(1, 0): 1}
Y: Poly = {(0, 1): 1}


def chi(n: int, var: Poly) -> Poly:
    """chi_0 = 1, chi_1 = var, chi_(n+1) = var*chi_n - chi_(n-1)."""
    a, b = ONE, var
    if n == 0:
        return a
    for _ in range(n - 1):
        a, b = b, padd(pmul(var, b), a, -1)
    return b


def S(b: int) -> Poly:
    return padd(chi(b, X), chi(b, Y))


def D(b: int) -> Poly:
    return padd(chi(b, X), chi(b, Y), -1)


def main() -> int:
    UMAX, AMAX = 24, 48
    failures: list[str] = []

    # 1. Clebsch-Gordan square decomposition.
    cg = 0
    for u in range(UMAX + 1):
        lhs = pmul(chi(u, X), chi(u, X))
        rhs: Poly = {}
        for i in range(u + 1):
            rhs = padd(rhs, chi(2 * i, X))
        if lhs != rhs:
            failures.append(f"chi_{u}^2 != sum chi_2i")
        else:
            cg += 1
    print(f"CLEBSCH_GORDAN_SQUARES checked={cg} max_u={UMAX}")

    # 2. The one known instance.
    if D(2) != pmul(D(1), S(1)):
        failures.append("D_2 != D_1 S_1")
    else:
        print("KNOWN_INSTANCE D_2 = D_1 S_1 PASS")

    # 3. Uniqueness over the range.
    tested = 0
    extra: list[str] = []
    for a in range(1, AMAX + 1):
        Da = D(a)
        for u in range(1, AMAX + 1):
            if (a, u) == (2, 1):
                continue
            tested += 1
            if Da == pmul(D(u), S(u)):
                extra.append(f"D_{a} = D_{u} S_{u}")
    if extra:
        failures.extend(extra)
    print(f"UNIQUENESS pairs_tested={tested} max_a={AMAX} extra_instances={len(extra)}")

    if failures:
        for f in failures:
            print(f"FAILURE: {f}", file=sys.stderr)
        print("SU2_DIFFERENCE_FACTORIZATION_LIMIT: FAILED", file=sys.stderr)
        return 1

    print(
        "SU2_DIFFERENCE_FACTORIZATION_LIMIT PASS "
        f"clebsch_gordan={cg} pairs={tested} instances=1 unique_pair=(a=2,u=1)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
