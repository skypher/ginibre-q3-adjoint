#!/usr/bin/env python3
"""Rational rounding and interval verification of AM-GM allocations.

The proposal stage (su2_amgm_lp) is floating point.  This module converts a
proposal into a rational allocation with a fixed denominator and then checks
every certificate condition again in mpmath interval arithmetic at 512 bits,
rebuilding the spectral data from scratch in interval mode.

The rebuild matters: verification must not inherit the proposal's floating
point values, or it would only be checking its own arithmetic.  Node values
are enclosed as 2cos(pi t/n) with interval pi and interval cosine, so every
enclosure is rigorous.

A certificate is accepted only when the interval lower bound of each slack is
strictly positive, or when the two sides are provably equal by structure.
"""

from __future__ import annotations

from fractions import Fraction

import mpmath as mp

VERIFY_PREC = 512


def _iv():
    ctx = mp.iv
    ctx.prec = VERIFY_PREC
    return ctx


def iv_nodes_and_characters(rank: int):
    """Rigorous interval enclosures of the weights and character table."""
    ctx = _iv()
    n = 2 * rank + 1
    labels = rank - 1
    weights, table = [], []
    for t in range(1, rank + 1):
        a = ctx.pi * t / n
        s = ctx.sin(a)
        weights.append(4 * s ** 2 / n)
        table.append([ctx.sin((2 * j + 1) * a) / s for j in range(1, labels + 1)])
    return weights, table


def iv_terms(signs, powers, weights, table):
    """Interval version of the chamber term list."""
    ctx = _iv()
    rank = len(weights)
    labels = len(signs)
    active = [l for l in range(labels) if signs[l] != 0]
    out = []
    for i in range(rank):
        for j in range(i + 1, rank):
            sign = 1
            coeff = 2 * weights[i] * weights[j]
            lam = [ctx.mpf(1)] * labels
            bad = False
            for l in active:
                b = table[i][l] + signs[l] * table[j][l]
                # A term is kept only when its base is bounded away from zero,
                # decided by the interval, not by a floating point tolerance.
                if b.a <= 0 <= b.b:
                    bad = True
                    break
                if b.b < 0 and powers[l] % 2 == 1:
                    sign = -sign
                ab = -b if b.b < 0 else b
                coeff = coeff * ab ** powers[l]
                lam[l] = b * b
            if not bad:
                out.append((sign, coeff, lam))
    return out


def round_allocation(alpha, denom: int):
    """Round a float allocation to multiples of 1/denom with rows summing to 1."""
    P = len(alpha)
    X = len(alpha[0]) if P else 0
    out = [[Fraction(0) for _ in range(X)] for _ in range(P)]
    for xi in range(X):
        units = []
        for pi in range(P):
            units.append(int(round(alpha[pi][xi] * denom)))
        total = sum(units)
        # Force the column to sum to exactly `denom` units.
        if total != denom:
            order = sorted(range(P), key=lambda pi: -alpha[pi][xi])
            k = 0
            while total != denom and k < 10 * P + 10:
                pi = order[k % P]
                step = 1 if total < denom else -1
                if units[pi] + step >= 0:
                    units[pi] += step
                    total += step
                k += 1
            if total != denom:
                return None
        for pi in range(P):
            out[pi][xi] = Fraction(units[pi], denom)
    return out


def verify(signs_chamber, powers, free, rank, alloc, pos, neg, denom):
    """Re-verify a rational allocation from scratch in interval arithmetic.

    Returns (ok, min_slack_lower_bound).
    """
    ctx = _iv()
    weights, table = iv_nodes_and_characters(rank)
    terms = iv_terms(signs_chamber, powers, weights, table)
    if len(terms) <= max(max(pos, default=-1), max(neg, default=-1)):
        return False, None

    # Absorb the regime's minimum exponent, exactly as the real problem does.
    lam = []
    coeff = []
    for (s, c, lm) in terms:
        cc = c
        for l in free:
            cc = cc * lm[l]
        coeff.append(cc)
        lam.append([lm[l] for l in free])

    # Rational structural conditions, checked exactly.
    for xi in range(len(neg)):
        if sum(alloc[pi][xi] for pi in range(len(pos))) != 1:
            return False, None
    for pi in range(len(pos)):
        if sum(alloc[pi][xi] for xi in range(len(neg))) > 1:
            return False, None

    log_lam = [[ctx.log(v) for v in row] for row in lam]
    log_c = [ctx.log(v) for v in coeff]

    worst = None
    for xi, x in enumerate(neg):
        # (D) domination in each free coordinate.
        for l in range(len(free)):
            acc = ctx.mpf(0)
            for pi, p in enumerate(pos):
                a = alloc[pi][xi]
                if a:
                    acc = acc + log_lam[p][l] * ctx.mpf(a.numerator) / a.denominator
            slack = acc - log_lam[x][l]
            if not (slack.a > 0 or slack.a == 0 == slack.b):
                return False, None
            worst = slack.a if worst is None else min(worst, slack.a)
        # (K) the constant condition.
        acc = ctx.mpf(0)
        for pi, p in enumerate(pos):
            a = alloc[pi][xi]
            if a:
                acc = acc + log_c[p] * ctx.mpf(a.numerator) / a.denominator
        slack = acc - log_c[x]
        if not (slack.a > 0 or slack.a == 0 == slack.b):
            return False, None
        worst = slack.a if worst is None else min(worst, slack.a)
    return True, worst
