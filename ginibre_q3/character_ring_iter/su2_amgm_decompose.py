#!/usr/bin/env python3
"""Recursive orthant decomposition for regimes no single AM-GM allocation closes.

At rank six, 408 residual regimes are genuinely LP-infeasible: no capacitated
allocation exists, so the regime has to be split before it can be certified.
This is the decomposition stage the original O11 proof referred to as
recursive orthant decompositions and translated orthants with exact faces.

The split
---------
Represent a node as terms `(s_t, c_t, lambda_t)` over free coordinates `F`,
meaning

    sum_t s_t c_t prod_(l in F) lambda_(t,l)^(u_l),      u_l >= 0.

For any `l in F` the index set splits exactly:

    {u_l >= 0} = {u_l = 0}  disjoint-union  {u_l >= 1}.

  * face `u_l = 0`   : drop `l` from `F`, coefficients unchanged;
  * tail `u_l >= 1`  : substitute `u_l = 1 + u'`, so `c_t *= lambda_(t,l)`
                       and `F` is unchanged.

Both branches are strictly simpler than the parent -- the face loses a
coordinate, and the tail strictly increases every coefficient's exponent -- and
their union is the parent exactly, so the decomposition is sound by
construction rather than by estimate.

A leaf with no free coordinates is a single exponent point, evaluated in
interval arithmetic.

Nothing here weakens the acceptance test: every node is closed either by an
interval-verified AM-GM certificate or by an interval evaluation of the actual
signed sum.
"""

from __future__ import annotations

import mpmath as mp

from su2_amgm_certify import _iv, iv_nodes_and_characters, iv_terms
from su2_amgm_lp import propose_allocation
from su2_amgm_certify import round_allocation


def _ladder(base: int):
    return [base, 2 * base, 5 * base, 10 * base, 50 * base]


def _leaf_nonnegative(iv_sign, iv_coeff) -> bool:
    """Evaluate a zero-dimensional leaf: is sum_t s_t c_t >= 0?"""
    ctx = _iv()
    total = ctx.mpf(0)
    for s, c in zip(iv_sign, iv_coeff):
        total = total + c if s > 0 else total - c
    return total.a >= 0


def _amgm_node(fl_sign, fl_coeff, fl_lam, iv_sign, iv_coeff, iv_lam, denom):
    """Try a single AM-GM allocation at one node.

    Floating point proposes from `fl_*`; the decision is taken on `iv_*` in
    interval arithmetic, exactly as in the flat stack.
    """
    ctx = _iv()
    nfree = len(fl_lam[0]) if fl_lam else 0
    if not any(s < 0 for s in fl_sign):
        return True                      # no negative term: nothing to dominate
    log_lam = [[mp.log(v) for v in row] for row in fl_lam]
    log_c = [mp.log(v) for v in fl_coeff]
    delta, proposal = propose_allocation(log_lam, log_c, fl_sign, nfree)
    if proposal is None or delta is None or delta <= 0:
        return False
    pos, neg, alpha = proposal

    iv_log_lam = [[ctx.log(v) for v in row] for row in iv_lam]
    iv_log_c = [ctx.log(v) for v in iv_coeff]

    for d in _ladder(denom):
        alloc = round_allocation(alpha, d)
        if alloc is None:
            continue
        ok = True
        for xi in range(len(neg)):
            if sum(alloc[pi][xi] for pi in range(len(pos))) != 1:
                ok = False
                break
        if ok:
            for pi in range(len(pos)):
                if sum(alloc[pi][xi] for xi in range(len(neg))) > 1:
                    ok = False
                    break
        if not ok:
            continue
        for xi, x in enumerate(neg):
            for l in range(nfree):
                acc = ctx.mpf(0)
                for pi, p in enumerate(pos):
                    a = alloc[pi][xi]
                    if a:
                        acc = acc + iv_log_lam[p][l] * ctx.mpf(a.numerator) / a.denominator
                sl = acc - iv_log_lam[x][l]
                if not (sl.a > 0 or sl.a == 0 == sl.b):
                    ok = False
                    break
            if not ok:
                break
            acc = ctx.mpf(0)
            for pi, p in enumerate(pos):
                a = alloc[pi][xi]
                if a:
                    acc = acc + iv_log_c[p] * ctx.mpf(a.numerator) / a.denominator
            sl = acc - iv_log_c[x]
            if not (sl.a > 0 or sl.a == 0 == sl.b):
                ok = False
                break
        if ok:
            return True
    return False


def certify(fl_sign, fl_coeff, fl_lam, iv_sign, iv_coeff, iv_lam,
            denom=100, depth=0, max_depth=6, budget=None):
    """Close a node by AM-GM, or split and recurse.  Returns (ok, stats)."""
    if budget is None:
        budget = {"nodes": 0, "leaves": 0, "amgm": 0, "cap": 4000}
    budget["nodes"] += 1
    if budget["nodes"] > budget["cap"]:
        return False, budget

    nfree = len(fl_lam[0]) if fl_lam else 0

    if nfree == 0:
        budget["leaves"] += 1
        return _leaf_nonnegative(iv_sign, iv_coeff), budget

    if _amgm_node(fl_sign, fl_coeff, fl_lam, iv_sign, iv_coeff, iv_lam, denom):
        budget["amgm"] += 1
        return True, budget

    if depth >= max_depth:
        return False, budget

    # Split on the coordinate where the negative terms hold the largest
    # advantage, since that is the obstruction the allocation could not pay.
    best_l, best_score = 0, None
    for l in range(nfree):
        worst = None
        for t in range(len(fl_sign)):
            if fl_sign[t] < 0:
                v = fl_lam[t][l]
                worst = v if worst is None else max(worst, v)
        score = worst if worst is not None else mp.mpf(0)
        if best_score is None or score > best_score:
            best_l, best_score = l, score

    l = best_l
    # Face u_l = 0: drop the coordinate.
    f_sign = list(fl_sign)
    f_coeff = list(fl_coeff)
    f_lam = [[row[j] for j in range(nfree) if j != l] for row in fl_lam]
    i_lam = [[row[j] for j in range(nfree) if j != l] for row in iv_lam]
    ok, budget = certify(f_sign, f_coeff, f_lam, iv_sign, iv_coeff, i_lam,
                         denom, depth + 1, max_depth, budget)
    if not ok:
        return False, budget

    # Tail u_l >= 1: absorb one factor of lambda into every coefficient.
    t_coeff = [fl_coeff[t] * fl_lam[t][l] for t in range(len(fl_coeff))]
    it_coeff = [iv_coeff[t] * iv_lam[t][l] for t in range(len(iv_coeff))]
    return certify(fl_sign, t_coeff, fl_lam, iv_sign, it_coeff, iv_lam,
                   denom, depth + 1, max_depth, budget)


def build_node(signs, powers, free, rank, weights, table):
    """Float and interval term data for a regime, ready for `certify`."""
    from su2_orbit_amgm import chamber_terms
    fl = chamber_terms(signs, powers, weights, table)
    fl_sign = [t.sign for t in fl]
    fl_coeff = []
    fl_lam = []
    for t in fl:
        c = t.coeff
        for l in free:
            c = c * t.lam[l]
        fl_coeff.append(c)
        fl_lam.append([t.lam[l] for l in free])

    ivw, ivt = iv_nodes_and_characters(rank)
    iv = iv_terms(signs, powers, ivw, ivt)
    iv_sign = [s for (s, _c, _l) in iv]
    iv_coeff = []
    iv_lam = []
    for (_s, c, lm) in iv:
        cc = c
        for l in free:
            cc = cc * lm[l]
        iv_coeff.append(cc)
        iv_lam.append([lm[l] for l in free])
    return fl_sign, fl_coeff, fl_lam, iv_sign, iv_coeff, iv_lam
