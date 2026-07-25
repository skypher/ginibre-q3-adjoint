#!/usr/bin/env python3
"""Rank-parameterised AM-GM certificate stack for odd orbit rings.

For a given rank it enumerates every support/sign/parity chamber and every
residual regime, skips the ones a direct capacitated Hall transport already
settles, and attempts an AM-GM certificate on the rest.

Pipeline per regime:
  1. propose a real allocation by linear programming (floating point);
  2. round it to a fixed rational denominator;
  3. re-verify from scratch in 512-bit interval arithmetic.

Only stage 3 decides.  A regime is reported closed only if its rational
allocation verifies.

Usage:
  verify_su2_orbit_amgm_stack.py --rank 6 [--denominator 100] [--limit N]
"""

from __future__ import annotations

import argparse
import sys
import time

import mpmath as mp

from su2_orbit_amgm import chamber_terms, nodes_and_characters
from su2_amgm_lp import propose_allocation
from su2_amgm_certify import round_allocation, verify

TOL = mp.mpf(10) ** -12


def denom_ladder(base: int):
    """Denominators to try, coarsest first.

    The recorded O11 run used denominator 100 throughout and 200 once, so the
    ladder starts there and escalates only when a margin is genuinely tight.
    """
    seen, out = set(), []
    for d in (base, 2 * base, 5 * base, 10 * base, 50 * base):
        if d not in seen:
            seen.add(d)
            out.append(d)
    return out


def hall_ok(terms, free) -> bool:
    """Direct capacitated Hall transport, matching the C++ census."""
    neg = [i for i, t in enumerate(terms) if t.sign < 0]
    if not neg:
        return True
    pos = [i for i, t in enumerate(terms) if t.sign > 0]
    mag = [t.coeff * mp.fprod([t.lam[l] for l in free]) if free else t.coeff
           for t in terms]

    def ge(a, b):
        return a + TOL * (1 + abs(a) + abs(b)) >= b

    edge = []
    for x in neg:
        bits = 0
        for k, p in enumerate(pos):
            if all(ge(terms[p].lam[l], terms[x].lam[l]) for l in free):
                bits |= 1 << k
        edge.append(bits)
    for s in range(1, 1 << len(neg)):
        d = mp.mpf(0)
        nb = 0
        for i in range(len(neg)):
            if (s >> i) & 1:
                d += mag[neg[i]]
                nb |= edge[i]
        cap = sum((mag[pos[k]] for k in range(len(pos)) if (nb >> k) & 1), mp.mpf(0))
        if not ge(cap, d):
            return False
    return True


def regimes(rank):
    """Yield (signs, powers, active, free) for every residual regime."""
    labels = rank - 1
    weights, table = nodes_and_characters(rank)
    for support in range(3 ** labels):
        z, signs = support, []
        for _ in range(labels):
            d = z % 3
            z //= 3
            signs.append(0 if d == 0 else (1 if d == 1 else -1))
        active = [l for l in range(labels) if signs[l]]
        if not active or not any(s < 0 for s in signs):
            continue
        for parity in range(1 << len(active)):
            powers = [0] * labels
            for i, l in enumerate(active):
                powers[l] = 1 if (parity >> i) & 1 else 2
            if sum(1 for l in active if signs[l] < 0 and powers[l] == 1) % 2:
                continue
            terms = chamber_terms(signs, powers, weights, table)
            if not terms or not any(t.sign < 0 for t in terms):
                continue
            for residual in range(1 << len(active)):
                free = [active[i] for i in range(len(active)) if (residual >> i) & 1]
                yield support, parity, signs, powers, terms, free


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", type=int, default=6)
    ap.add_argument("--denominator", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0, help="stop after N residual regimes")
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    rank = args.rank
    denom = args.denominator
    t0 = time.time()

    direct = closed = failed_lp = failed_round = failed_verify = 0
    denom_used: dict[int, int] = {}
    seen = 0
    worst_margin = None
    unresolved = []

    for support, parity, signs, powers, terms, free in regimes(rank):
        if hall_ok(terms, free):
            direct += 1
            continue
        seen += 1
        if args.limit and seen > args.limit:
            seen -= 1
            break

        log_lam = [[mp.log(t.lam[l]) for l in free] for t in terms]
        log_c = [mp.log(t.coeff) for t in terms]
        tsigns = [t.sign for t in terms]

        delta, proposal = propose_allocation(log_lam, log_c, tsigns, len(free))
        if proposal is None or delta is None or delta <= 0:
            failed_lp += 1
            unresolved.append((support, parity, tuple(free), "lp"))
            continue
        pos, neg, alpha = proposal
        # A coarse denominator can eat the LP margin.  Escalate until the
        # rounded allocation verifies, and record the denominator actually
        # needed so the certificate states its own cost.
        ok, margin, used = False, None, None
        for d in denom_ladder(denom):
            alloc = round_allocation(alpha, d)
            if alloc is None:
                continue
            ok, margin = verify(signs, powers, free, rank, alloc, pos, neg, d)
            if ok:
                used = d
                break
        if not ok:
            failed_verify += 1
            unresolved.append((support, parity, tuple(free), "verify"))
            continue
        closed += 1
        denom_used[used] = denom_used.get(used, 0) + 1
        if margin is not None:
            worst_margin = margin if worst_margin is None else min(worst_margin, margin)
        if args.progress and closed % 100 == 0:
            print(f"  ... {closed} closed / {seen} attempted", file=sys.stderr)

    elapsed = time.time() - t0
    level = 2 * rank - 1
    print(f"SU2_ORBIT_AMGM_STACK rank={rank} level={level} denominator={denom}")
    print(f"  direct_hall      {direct}")
    print(f"  residual         {seen}")
    print(f"  closed           {closed}")
    print(f"  lp_infeasible    {failed_lp}")
    print(f"  round_failed     {failed_round}")
    if denom_used:
        blurb = " ".join(f"{d}:{c}" for d, c in sorted(denom_used.items()))
        print(f"  denominators     {blurb}")
    print(f"  verify_failed    {failed_verify}")
    if worst_margin is not None:
        print(f"  min_slack        {mp.nstr(mp.mpf(worst_margin), 8)}")
    print(f"  elapsed_seconds  {elapsed:.1f}")
    if seen:
        print(f"  coverage         {100.0*closed/seen:.1f}% of residual regimes")
    for u in unresolved[:20]:
        print(f"UNRESOLVED support={u[0]} parity={u[1]} free={u[2]} stage={u[3]}")
    if len(unresolved) > 20:
        print(f"UNRESOLVED ... and {len(unresolved)-20} more")
    return 0


if __name__ == "__main__":
    sys.exit(main())
