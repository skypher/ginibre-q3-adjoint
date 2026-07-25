#!/usr/bin/env python3
"""Rank-parameterised AM-GM certificate stack for the odd simple-current
orbit rings `O_(2m-1)`.

Background
----------
Each orbit ring's residual problem is: for a support/sign/parity chamber and
a choice of free coordinates, show

    sum_t s_t c_t prod_l lambda_(t,l)^(r_l) >= 0        for all r_l >= 0,

where `s_t` is a sign, `c_t>0`, and `lambda_(t,l)>0`.  Terms come from
off-diagonal spectral node pairs; diagonal pairs vanish whenever the word has
a minus factor.

Certificate
-----------
Allocate a share `alpha[p][x]` of each positive term `p` to each negative term
`x`.  Applying weighted AM-GM with the allocation itself as the weight vector
turns the requirement into four families of *linear* conditions on alpha:

    (N)  sum_p alpha[p][x] = 1                         for each negative x
    (C)  sum_x alpha[p][x] <= 1                        for each positive p
    (D)  sum_p alpha[p][x] log lambda[p][l]
             >= log lambda[x][l]                       for each x, free l
    (K)  sum_p alpha[p][x] log c[p] >= log c[x]        for each negative x

So certificate search is a linear feasibility problem -- a capacitated
transportation problem -- not a general nonlinear search.  Floating point is
used only to *propose* an allocation; the proposal is rounded to a rational
with fixed denominator and then re-checked in interval arithmetic with
directed rounding.  No floating-point decision is accepted as evidence.

This module is the shared engine.  It is parameterised by rank, so the same
code serves `O_9`, `O_11`, `O_13` and beyond.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable

import mpmath as mp

# Working precision for the proposal stage.
mp.mp.prec = 240
# Independent, higher precision for the interval verification stage.
VERIFY_PREC = 512


class Term:
    """A single spectral node-pair term."""

    __slots__ = ("sign", "coeff", "lam")

    def __init__(self, sign: int, coeff, lam: list):
        self.sign = sign
        self.coeff = coeff
        self.lam = lam

    def __repr__(self) -> str:
        return f"Term(sign={self.sign:+d}, coeff={mp.nstr(self.coeff,6)})"


def nodes_and_characters(rank: int):
    """Return (weights, character_table) for the orbit ring of the given rank.

    Level is k=2*rank-1 and the Verlinde order is n=k+2=2*rank+1.  Node t
    carries x_t=2cos(pi t/n); the orbit character B_j equals chi_(2j).  This
    matches analyze_su2_o13_frontier.cpp, which was validated against the
    recorded O11 census.
    """
    n = 2 * rank + 1
    labels = rank - 1
    weights, table = [], []
    for t in range(1, rank + 1):
        a = mp.pi * t / n
        weights.append(4 * mp.sin(a) ** 2 / n)
        row = [mp.sin((2 * j + 1) * a) / mp.sin(a) for j in range(1, labels + 1)]
        table.append(row)
    return weights, table


def chamber_terms(signs: list[int], powers: list[int], weights, table) -> list[Term]:
    """Build the off-diagonal term list for one support/parity chamber.

    `signs[l]` is 0, +1 or -1 for label B_(l+1); `powers[l]` is that label's
    minimum exponent (1 or 2).  Returns [] when the chamber is degenerate.
    """
    rank = len(weights)
    labels = len(signs)
    active = [l for l in range(labels) if signs[l] != 0]
    terms: list[Term] = []
    for i in range(rank):
        for j in range(i + 1, rank):
            sign = 1
            coeff = 2 * weights[i] * weights[j]
            lam = [mp.mpf(1)] * labels
            degenerate = False
            for l in active:
                b = table[i][l] + signs[l] * table[j][l]
                if abs(b) < mp.mpf(10) ** (-40):
                    degenerate = True
                    break
                q = powers[l]
                if b < 0 and (q % 2 == 1):
                    sign = -sign
                coeff *= abs(b) ** q
                lam[l] = b * b
            if not degenerate:
                terms.append(Term(sign, coeff, lam))
    return terms


def restrict(terms: list[Term], free: Iterable[int]) -> tuple[list[Term], list[int]]:
    """Fix non-free coordinates at their minimum and keep free ones.

    In a regime the free coordinates satisfy r_l>=1, so absorb one factor of
    lambda into the coefficient and let the remaining exponent run from zero.
    """
    free = list(free)
    out = []
    for t in terms:
        c = t.coeff
        for l in free:
            c = c * t.lam[l]
        out.append(Term(t.sign, c, [t.lam[l] for l in free]))
    return out, free
