#!/usr/bin/env python3
"""Order-cone fan decomposition for regimes with a diagonal obstruction.

Axis-aligned orthant splitting cannot help a regime whose obstruction lies
along a diagonal: a negative term that is dominated when one coordinate grows
alone, but not when several grow together.  Those regimes stay LP-infeasible
however far the face/tail recursion goes.

The fan
-------
For free coordinates `u_1..u_k >= 0`, the `k!` order cones

    C_sigma = { u_(sigma(1)) >= u_(sigma(2)) >= ... >= u_(sigma(k)) >= 0 }

cover the orthant.  Each is unimodular, with the exact integer substitution

    u_(sigma(j)) = sum_(i >= j) m_i,          m_i >= 0,

whose inverse is the difference of consecutive coordinates, so lattice points
of the cone correspond exactly to nonnegative integer `m`.  Substituting,

    prod_l lambda_(t,l)^(u_l) = prod_i ( prod_(j <= i) lambda_(t,sigma(j)) )^(m_i),

so the cone is again an AM-GM problem in `m`, with bases replaced by running
products along the order.  Coefficients are untouched because every cone has
its apex at the origin.

That means the same linear program certifies each cone; only the base vectors
change.  A regime is closed when every one of the `k!` cones is closed.

Covering, not partitioning, is all that is needed: cones meet on their
boundaries and a lattice point lying in several is simply certified more than
once.
"""

from __future__ import annotations

from itertools import permutations

from su2_amgm_decompose import _amgm_node


def running_products(lam, order):
    """Bases for one order cone: prod_(j<=i) lambda[sigma(j)]."""
    out = []
    for row in lam:
        acc, vals = None, []
        for i in range(len(order)):
            acc = row[order[i]] if acc is None else acc * row[order[i]]
            vals.append(acc)
        out.append(vals)
    return out


def certify_fan(fl_sign, fl_coeff, fl_lam, iv_sign, iv_coeff, iv_lam,
                denom=100, max_cones=720):
    """Close a regime by certifying every order cone.  Returns (ok, cones)."""
    k = len(fl_lam[0]) if fl_lam else 0
    if k <= 1:
        # One coordinate: the fan is the orthant itself, nothing gained.
        return False, 0
    orders = list(permutations(range(k)))
    if len(orders) > max_cones:
        return False, 0
    done = 0
    for order in orders:
        f_lam = running_products(fl_lam, order)
        i_lam = running_products(iv_lam, order)
        if not _amgm_node(fl_sign, fl_coeff, f_lam,
                          iv_sign, iv_coeff, i_lam, denom):
            return False, done
        done += 1
    return True, done
