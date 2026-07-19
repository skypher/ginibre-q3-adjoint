#!/usr/bin/env python3
"""Fit Laplace-expansion coefficients c_1, c_2, c_3 from our exact m_n data.

Model: m_n = tilde m_n * (r_infty + c_1/n + c_2/n^2 + c_3/n^3 + ...)
where tilde m_n = d^n * n^{-d/2} * exp(-d^2/(48n)).

The fit uses least squares over our top 8-10 values of n for each family.
Non-zero r_infty absorbs the BPV prefactor K_G.

Results per family (from runs on optimus):
  G2 (d=14): c_1 ≈ -12, c_2 ≈ 88, c_3 ≈ -480
  F4 (d=52): c_1 ≈ -109, c_2 ≈ 5009, c_3 ≈ -92700
  E6 (d=78): c_1 ≈ -151, c_2 ≈ 8267, c_3 ≈ -159420
  E7 (d=133): c_1 ≈ -171, c_2 ≈ 9872, c_3 ≈ -192145
  E8 (d=248): c_1 ≈ -186, c_2 ≈ 11524, c_3 ≈ -238613

Observed pattern: |c_k|/|c_1|^k is decreasing, but sub-asymptotic at our
verified ranges (c_2/n^2 is comparable to c_1/n for n much smaller than
|c_1|).
"""
import sys, os, re, math
import numpy as np

def parse(path):
    moms = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'm_(\d+) = (\d+)', line)
            if m:
                moms[int(m.group(1))] = int(m.group(2))
    return moms

def fit_coefficients(log_path, d, order=3):
    moms = parse(log_path)
    K = max(moms.keys())
    # Use last ~10 values for fit
    n_start = max(K - 9, 10)
    n_list = list(range(n_start, K + 1))
    # log tilde_m = n log d - (d/2) log n - d^2/(48 n)
    log_y = []
    for n in n_list:
        lm = math.log(moms[n])
        lt = n * math.log(d) - (d/2) * math.log(n) - d*d/(48*n)
        log_y.append(lm - lt)
    y = np.array([math.exp(ly - log_y[-1]) for ly in log_y])
    inv = np.array([1.0/n for n in n_list])
    X = np.column_stack([inv**j for j in range(order + 1)])
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    r_inf = coeffs[0]
    cs = [coeffs[j+1] / r_inf for j in range(order)]
    return cs, r_inf

if __name__ == '__main__':
    groups = {
        'G2': (14, 'logs/g2_200.log'),
        'F4': (52, 'logs/f4_100.log'),
        'E6': (78, 'logs/e6_80.log'),
        'E7': (133, 'logs/e7_70.log'),
        'E8': (248, 'logs/e8_70.log'),
    }
    here = os.path.dirname(os.path.abspath(__file__))
    print(f"{'G':<4}{'d':>5}{'c_1':>12}{'c_2':>12}{'c_3':>14}")
    print('-' * 50)
    for G, (d, log) in groups.items():
        path = os.path.join(here, log)
        cs, _ = fit_coefficients(path, d, order=3)
        print(f"{G:<4}{d:>5}{cs[0]:>12.1f}{cs[1]:>12.1f}{cs[2]:>14.1f}")
