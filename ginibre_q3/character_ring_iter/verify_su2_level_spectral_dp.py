#!/usr/bin/env python3
"""Cross-validate the full-level spectral corner against exact fusion DP.

The certificate engine is being extended from odd orbit rings to the full
fusion ring SU(2)_k at arbitrary level, using

    J = sum_{r,s} w_r w_s prod_i ( chi_{a_i}(r) + eps_i chi_{a_i}(s) ),

with the Verlinde normalisation of the level 1-4 note:

    w_j    = 2/(k+2) sin^2((j+1)pi/(k+2)),
    chi_a(j) = sin((a+1)(j+1)pi/(k+2)) / sin((j+1)pi/(k+2)).

This script checks that formula end to end against exact integer dynamic
programming in the doubled fusion ring: J is the [V_0 tensor V_0] coefficient
of prod_i (V_{a_i} tensor V_0 + eps_i V_0 tensor V_{a_i}), an integer.

Words are drawn at random with no support-disjointness restriction and with
both parities of the total degree, so the check also confirms the parity
selection rule (odd total degree gives exactly zero) on the spectral side.
"""

from __future__ import annotations

import math
import random
import sys


def fusion_range(a: int, b: int, k: int):
    return range(abs(a - b), min(a + b, 2 * k - a - b) + 1, 2)


def exact_corner(word, k: int) -> int:
    """Integer J by DP over coefficient matrices in the doubled ring."""
    dim = k + 1
    coeff = [[0] * dim for _ in range(dim)]
    coeff[0][0] = 1
    for (a, eps) in word:
        new = [[0] * dim for _ in range(dim)]
        for c in range(dim):
            for d in range(dim):
                v = coeff[c][d]
                if not v:
                    continue
                for c2 in fusion_range(a, c, k):
                    new[c2][d] += v
                for d2 in fusion_range(a, d, k):
                    new[c][d2] += eps * v
        coeff = new
    return coeff[0][0]


def spectral_corner(word, k: int) -> float:
    n = k + 2
    nodes = range(k + 1)
    w = [2.0 / n * math.sin((j + 1) * math.pi / n) ** 2 for j in nodes]
    chi = [[math.sin((a + 1) * (j + 1) * math.pi / n)
            / math.sin((j + 1) * math.pi / n)
            for a in range(k + 1)] for j in nodes]
    total = 0.0
    for r in nodes:
        for s in nodes:
            p = w[r] * w[s]
            for (a, eps) in word:
                p *= chi[r][a] + eps * chi[s][a]
            total += p
    return total


def main() -> int:
    rng = random.Random(20260725)
    trials = 0
    worst = 0.0
    for _ in range(240):
        k = rng.choice([2, 3, 4, 5, 6, 7, 8])
        length = rng.randint(1, 7)
        word = [(rng.randint(1, k), rng.choice([1, -1])) for _ in range(length)]
        exact = exact_corner(word, k)
        approx = spectral_corner(word, k)
        err = abs(approx - exact) / max(1.0, abs(exact))
        worst = max(worst, err)
        trials += 1
        if err > 1e-9:
            print(f"MISMATCH k={k} word={word} exact={exact} spectral={approx}",
                  file=sys.stderr)
            return 1
    print(f"SU2_LEVEL_SPECTRAL_DP PASS trials={trials} worst_rel_err={worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
