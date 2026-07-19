#!/usr/bin/env python3
"""For F_4, the Wick computation showed <I_6(F_4)>_mu = 0 exactly.

Question: does this come from a polynomial identity analogous to G_2's
  prod_{alpha in R_+(G_2)}(alpha.w)^2 = (27/1024)(|w|^12 - I_6^(0)(w)^2)?

Checks: (a) angular-avg of (I_6 candidate) on the sphere S^3 is zero, OR
        (b) some other identity forces vanishing.

Candidate primary I_6(F_4) := bar_R - mu_1 |w|^6 where
  mu_1 = B_3 * 15 / (r(r+2)(r+4)) = 108 * 15 / (4*6*8) = 1620/192 = 135/16.

We numerically sample the sphere S^3 and check if
  int_{S^3} I_6_candidate(hat w) prod(alpha.hat w)^2 d Omega
vanishes.  If yes, that would explain <I_6>_mu = 0 structurally
(since the measure factorizes radial times angular).
"""
import numpy as np
from fractions import Fraction

def f4_positive_roots():
    roots = []
    for i in range(4):
        for j in range(i+1, 4):
            for s in (1, -1):
                v = [0.0] * 4
                v[i] = 1.0
                v[j] = float(s)
                roots.append(v)
    for i in range(4):
        v = [0.0] * 4
        v[i] = 1.0
        roots.append(v)
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                roots.append([0.5, s1/2, s2/2, s3/2])
    return np.array(roots, dtype=np.float64)

def main():
    roots = f4_positive_roots()
    assert len(roots) == 24
    r, Rp, d = 4, 24, 52
    mu_1 = 108 * 15 / (r * (r+2) * (r+4))  # 135/16 = 8.4375
    print(f"F_4: mu_1 = {mu_1}")

    # Sample uniformly on S^3: draw from N(0,I_4) and normalize
    rng = np.random.default_rng(0)
    N = 5_000_000
    total_I6_weighted = 0.0
    total_weight = 0.0
    batch = 100_000
    for _ in range(N // batch):
        W = rng.standard_normal(size=(batch, 4))
        norm = np.linalg.norm(W, axis=1, keepdims=True)
        hatW = W / norm
        # prod(alpha . hat w)^2
        inner = hatW @ roots.T
        weight_angular = np.prod(inner * inner, axis=1)
        # bar R = sum (alpha.hat w)^6
        barR = np.sum(inner**6, axis=1)
        # mu_1 |hat w|^6 = mu_1 (since unit vector)
        I6_cand = barR - mu_1  # evaluated on unit sphere
        total_I6_weighted += (I6_cand * weight_angular).sum()
        total_weight += weight_angular.sum()
    mean_I6_ang = total_I6_weighted / total_weight
    print(f"F_4: angular avg of I_6_cand * prod(alpha.hat w)^2, normalized = {mean_I6_ang:.4e}")
    print(f"F_4: should be ~0 if polynomial identity analogous to G_2 holds")

    # Also check just the angular avg of I_6 (no weight):
    total_I6_plain = 0.0
    for _ in range(N // batch):
        W = rng.standard_normal(size=(batch, 4))
        norm = np.linalg.norm(W, axis=1, keepdims=True)
        hatW = W / norm
        inner = hatW @ roots.T
        barR = np.sum(inner**6, axis=1)
        I6 = barR - mu_1
        total_I6_plain += I6.sum()
    mean_I6_plain = total_I6_plain / N
    print(f"F_4: UNWEIGHTED angular avg of I_6_cand = {mean_I6_plain:.4e}")
    print(f"F_4: should be ~0 by construction (mu_1 is chosen so this vanishes)")

if __name__ == '__main__':
    main()
