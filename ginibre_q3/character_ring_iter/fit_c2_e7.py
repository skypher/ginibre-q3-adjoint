#!/usr/bin/env python3
"""Higher-order fit of c_2(E_7) from character-ring moments, with
c_1 fixed at closed form.  Fits c_2, c_3, c_4, c_5 simultaneously to
log(m_n / tilde_m_n) - c_1/n via least-squares.

The goal: cross-check the MC value c_2(E_7) ~ 680351 against a direct
fit of the character-ring moments that uses high-order expansion.
If fit and MC agree, both provide converging evidence.

The low-order fit c_2 ~ 563000 (paper's empirical) comes from truncating
at c_2/n^2 and fitting with n in [60, 70], which is deeply pre-asymptotic
for E_7 (|c_1| = 1139, so we need n >> 1139 for asymptotic-sum truncation).

Here we keep n in [10, 70] and fit up to c_5 to absorb pre-asymptotic
terms and extract c_2 robustly.
"""
import os, re, math
import numpy as np
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

def parse(path):
    moms = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'm_(\d+) = (\d+)', line)
            if m:
                moms[int(m.group(1))] = int(m.group(2))
    return moms

def fit_all_c(log_path, d, c1_closed, orders_max=5, n_min=10, n_max=None):
    moms = parse(log_path)
    K = max(moms.keys())
    if n_max is None:
        n_max = K
    ns = list(range(n_min, n_max + 1))
    ns = [n for n in ns if n in moms and moms[n] > 0]
    # log(m_n / tilde_m_n) where tilde_m_n = d^{n + d/2} n^{-d/2} * const
    log_m = [math.log(moms[n]) for n in ns]
    log_t = [n * math.log(d) + (d/2) * math.log(d) - (d/2) * math.log(n) for n in ns]
    # residual: log(m/tilde_m) = log(const) + c_1/n + c_2/n^2 + c_3/n^3 + ...
    # Note: log(1 + c_1/n + c_2/n^2 + ...) != c_1/n + c_2/n^2 exactly;
    # relating to tilde c's via
    # log(1 + x) = x - x^2/2 + ... where x = c_1/n + c_2/n^2 + ...
    # So the coefficients are tilde_c_k where
    #   tilde_c_1 = c_1
    #   tilde_c_2 = c_2 - c_1^2/2
    #   tilde_c_3 = c_3 - c_1 c_2 + c_1^3/3
    #   tilde_c_4 = c_4 - c_1 c_3 - c_2^2/2 + c_1^2 c_2 - c_1^4/4
    #   tilde_c_5 = (algebraic combination; increases in complexity)
    #
    # We fit the tilde_c_j's and back out c_2 from tilde_c_2.
    # Known: tilde_c_1 = c_1 (fixed).
    residual = [lm - lt - c1_closed / n for lm, lt, n in zip(log_m, log_t, ns)]
    # Build design matrix: columns [1, 1/n^2, 1/n^3, 1/n^4, 1/n^5, ...]
    # (1/n omitted since c_1 fixed; constant column for log_const; 1/n^2 onward)
    cols = [[1.0] * len(ns)] + [[1.0 / n**k for n in ns] for k in range(2, orders_max + 1)]
    X = np.array(cols).T
    y = np.array(residual, dtype=float)
    coeffs, resid_norm, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    log_const = coeffs[0]
    tilde_c = [c1_closed] + list(coeffs[1:])
    # Recover c_2 from tilde_c_2 = c_2 - c_1^2/2
    c_2 = tilde_c[1] + c1_closed**2 / 2
    return log_const, c_2, tilde_c, coeffs, ns

def sliding_window_fit(log_path, d, c1_closed, orders_max=5, window_size=20):
    """Try several n-ranges and report c_2 per window."""
    moms = parse(log_path)
    K = max(moms.keys())
    results = []
    for n_start in range(10, K - window_size + 1, 5):
        n_end = n_start + window_size - 1
        if n_end > K:
            continue
        try:
            lc, c2, tilde_c, coeffs, ns = fit_all_c(log_path, d, c1_closed,
                                                   orders_max=orders_max,
                                                   n_min=n_start, n_max=n_end)
            results.append((n_start, n_end, c2, tilde_c))
        except Exception as e:
            continue
    return results

def main():
    groups = {
        'G2': (14, 'logs/g2_200.log', 143.967882),  # exact c_2
        'F4': (52, 'logs/f4_200.log', 17287.4588),
        'E6': (78, 'logs/e6_80.log', 82866.74),  # MC value
        'E7': (133, 'logs/e7_70.log', 680350.58),  # MC value
        'E8': (248, 'logs/e8_70.log', 346377539/45),  # exact
    }
    print("Fit c_2 with c_1 fixed at closed form, higher-order terms free")
    for G, (d, log_rel, c2_ref) in groups.items():
        c1_closed = -d * (d + 4) / 16
        path = os.path.join(HERE, log_rel)
        print(f"\n=== {G} (d = {d}, c_1 = {c1_closed}, c_2_ref = {c2_ref}) ===")
        for orders in [3, 4, 5, 6]:
            try:
                lc, c2, tilde_c, coeffs, ns = fit_all_c(path, d, c1_closed,
                                                       orders_max=orders,
                                                       n_min=10)
                diff_pct = 100 * (c2 - c2_ref) / c2_ref if c2_ref else 0
                print(f"  orders = {orders}: n in [{min(ns)}, {max(ns)}]  "
                      f"c_2_fit = {c2:.4f}  (diff from ref: {diff_pct:+.3f}%)")
            except Exception as e:
                print(f"  orders = {orders}: FAILED {e}")

if __name__ == '__main__':
    main()
