#!/usr/bin/env python3
"""Empirical shift-2 log-convexity analysis for all m_3 = 1 exceptional families.

Tests:
  (a) m_n m_{n+4} >= m_{n+2}^2 directly, for all n in the verified range.
  (b) Laplace-series residual R_n := m_n / tilde m_n - 1 - c_1/n - c_2/n^2,
      where tilde m_n = d^{n+d/2} n^{-d/2} * I_0/|W|  (unknown normalizing
      constant absorbed by least-squares fit of r_infty).
      Report R_n * n^3 as an estimate of c_3(n); if this is bounded, the
      Wong remainder at order 2 is effectively controlled.
  (c) Minimum n at which the Laplace third-order series is guaranteed
      positive (i.e., 1 + c_1/n + c_2/n^2 - C_3/n^3 >= 1/2) given an
      empirically-estimated C_3.

All arithmetic in Python's arbitrary-precision ints + mpmath for log.
"""
import os, re, math, sys
from fractions import Fraction

try:
    import mpmath as mp
    mp.mp.dps = 60
except ImportError:
    print("mpmath required; install with pip install mpmath", file=sys.stderr)
    sys.exit(1)

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = {
    'G2': (14, os.path.join(HERE, 'logs', 'g2_200.log')),
    'F4': (52, os.path.join(HERE, 'logs', 'f4_200.log')),
    'E6': (78, os.path.join(HERE, 'logs', 'e6_80.log')),
    'E7': (133, os.path.join(HERE, 'logs', 'e7_70.log')),
    'E8': (248, os.path.join(HERE, 'logs', 'e8_70.log')),
}

def parse_log(path):
    moms = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'm_(\d+) = (\d+)', line)
            if m:
                moms[int(m.group(1))] = int(m.group(2))
    return moms

def shift2_check(moms, n_start=4, n_stop=None):
    """Return list of (n, m_n m_{n+4} - m_{n+2}^2, sign)."""
    K = max(moms.keys())
    if n_stop is None:
        n_stop = K - 4
    result = []
    for n in range(n_start, n_stop + 1):
        if n + 4 not in moms or n + 2 not in moms or n not in moms:
            continue
        diff = moms[n] * moms[n+4] - moms[n+2]**2
        result.append((n, diff))
    return result

def laplace_residual(moms, d, n_min=10):
    """R_n := log(m_n) - [n log d + (d/2) log d - (d/2) log n + c_1/n + const].

    Using closed form c_1 = -d(d+4)/16.
    Returns list of (n, log(m_n) - model_without_c2_c3, R_n * n^3).
    """
    c1 = -d*(d+4) / 16
    K = max(moms.keys())
    log_tilde_over_const = []  # log(tilde m_n) - log(I_0/|W|)
    log_m = []
    ns = []
    for n in range(n_min, K + 1):
        if n not in moms or moms[n] <= 0:
            continue
        lm = mp.log(moms[n])
        lt = n * mp.log(d) + (d/2) * mp.log(d) - (d/2) * mp.log(n)
        log_m.append(lm)
        log_tilde_over_const.append(lt)
        ns.append(n)
    # Fit log const by least squares using 1 + c_1/n model
    residuals_with_const = []
    for (n, lm, lt) in zip(ns, log_m, log_tilde_over_const):
        # log(m/tilde_m) = log(I_0/|W|) + log(1 + c_1/n + c_2/n^2 + ...)
        # For large n, log(1 + c_1/n) ≈ c_1/n
        # So: lm - lt ≈ log const + c_1/n + higher
        resid = lm - lt
        residuals_with_const.append(resid)
    return ns, residuals_with_const, c1

def fit_laplace_asymptotic(ns, resids, c1, orders=(1, 2, 3)):
    """Fit resids = log_const + c_1/n + beta_2/n^2 + beta_3/n^3.
    But we hold c_1 fixed to closed-form. Fit log_const, beta_2, beta_3.
    Using last 10-15 points for the fit."""
    # Last few points, skipping pre-asymptotic
    k = min(15, len(ns))
    ns_fit = ns[-k:]
    rs_fit = resids[-k:]
    # Build matrix with col 0 = 1 (log_const), col 1 = 1/n^2, col 2 = 1/n^3
    # Known c_1/n subtracted from rs_fit
    import numpy as np
    rs_adj = [float(r) - c1 / n for r, n in zip(rs_fit, ns_fit)]
    cols = [[1.0] * k] + [[1.0 / (n ** p) for n in ns_fit] for p in [2, 3]]
    X = np.array(cols).T
    y = np.array(rs_adj)
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    log_const, beta_2, beta_3 = coeffs
    # Beta_2 corresponds to log(1 + c_2/n^2) ≈ c_2/n^2 - c_1 c_2/(2 n^3) + ...
    # So log const term absorbs prefactor, beta_2 ≈ c_2 (tilde), beta_3 ≈ c_3 - c_1 c_2 / 2
    return log_const, beta_2, beta_3

def analyze(name, d, log_path):
    print(f"\n{'='*70}")
    print(f"Family: {name}, d = {d}, c_1(closed form) = -{d*(d+4)//16}")
    print('='*70)
    moms = parse_log(log_path)
    K = max(moms.keys())
    print(f"Moments available: m_0 through m_{K} (total {len(moms)} entries)")

    # (a) shift-2 log-convexity check
    print(f"\n(a) Checking m_n m_{{n+4}} >= m_{{n+2}}^2 for n = 4..{K-4}")
    sh2 = shift2_check(moms, n_start=4, n_stop=K-4)
    violations_even = [(n, diff) for n, diff in sh2 if n % 2 == 0 and diff < 0]
    violations_odd = [(n, diff) for n, diff in sh2 if n % 2 == 1 and diff < 0]
    print(f"   even n: {len(violations_even)} violations out of "
          f"{sum(1 for n,_ in sh2 if n%2==0)} checks")
    print(f"   odd  n: {len(violations_odd)} violations out of "
          f"{sum(1 for n,_ in sh2 if n%2==1)} checks")
    if violations_odd:
        print(f"   First odd violations: {violations_odd[:3]}")
    else:
        print(f"   [all odd n log-convex, no violations]")

    # Specific small-n cases
    print(f"\n   Log-convexity ratio log10(m_n m_{{n+4}}/m_{{n+2}}^2) per odd n (sample):")
    odd_sh2 = [(n, diff) for n, diff in sh2 if n % 2 == 1]
    for (n, diff) in odd_sh2[:5] + ([] if len(odd_sh2) < 10 else [odd_sh2[-3]] + [odd_sh2[-1]]):
        mn, mn2, mn4 = moms[n], moms[n+2], moms[n+4]
        # diff / mn2^2 ≈ small number
        if mn2 > 0:
            ratio = float(mp.mpf(diff) / mp.mpf(mn2)**2)
            lr = mp.log10(mp.mpf(mn) * mp.mpf(mn4) / mp.mpf(mn2)**2)
        else:
            ratio = float('nan'); lr = 0
        print(f"    n={n:4}: m_n m_{{n+4}}/m_{{n+2}}^2 = 1 + {ratio:+.3e} (log10 = {float(lr):+.4e})")

    # (b) Laplace residual
    print(f"\n(b) Laplace-residual diagnostics (model 1 + c_1/n + c_2/n^2 + ...)")
    ns, resids, c1 = laplace_residual(moms, d, n_min=10)
    # Fit
    try:
        lc, b2, b3 = fit_laplace_asymptotic(ns, resids, c1)
        print(f"   Fit (last 15 ns, c_1 fixed = {c1}):")
        print(f"     log const: {lc:.4f}")
        print(f"     beta_2    (~ c_2 - c_1^2/2): {b2:.4f}   (expected c_2_{{paper}} - c_1^2/2)")
        print(f"     beta_3    (~ c_3 - c_1 c_2 / 2 + c_1^3/3): {b3:.4f}")
        # Compute tail residuals R_n * n^3
        print(f"\n   Residual R_n := log(m/tilde) - (log const + c_1/n + beta_2/n^2):")
        print(f"   n^3 * R_n = c_3 + O(1/n)  (reported; fluctuation indicates model quality)")
        for j in range(-min(15, len(ns)), 0, 2):
            n = ns[j]
            r = resids[j] - (lc + c1/n + b2/n**2)
            rn3 = float(r * n**3)
            print(f"     n = {n:4}: n^3 R_n = {rn3:+12.3e}")
    except Exception as e:
        print(f"   Fit failed: {e}")

    # (c) Effective threshold
    print(f"\n(c) Effective Laplace threshold (series 1 + c_1/n + c_2/n^2 >= 1/2)")
    c1_val = -d * (d+4) / 16
    # Need an explicit c_2 value per family
    c2_vals = {'G2': 144.0, 'F4': 17287.5, 'E6': 82150.0, 'E7': 563000.0, 'E8': 7697279.0/45*45}
    # (actually E8 is exact: 346_377_539/45)
    from fractions import Fraction
    c2_exact = {'E8': Fraction(346_377_539, 45)}
    c2 = c2_vals.get(name)
    if c2 is not None:
        # Solve: 1 + c_1/n + c_2/n^2 = 1/2  =>  n^2/2 + c_1 n + c_2 = 0  =>  n^2 + 2 c_1 n + 2 c_2 = 0
        # n = -c_1 ± sqrt(c_1^2 - 2 c_2)
        disc = c1_val**2 - 2 * c2
        if disc < 0:
            print(f"    c_1^2 - 2 c_2 = {disc:.4f} < 0: series always > 1/2 for n > 0")
            print(f"    Threshold from (2nd-order) Laplace alone is effectively 0.")
        else:
            n_plus = -c1_val + math.sqrt(disc)
            print(f"    1+c_1/n+c_2/n^2 = 1/2 has roots at n ≈ {-c1_val - math.sqrt(disc):.2f} and {n_plus:.2f}")
            print(f"    So 2nd-order series >= 1/2 for n >= {n_plus:.1f}, to compare with 2|c_1| = {-2*c1_val}")

if __name__ == '__main__':
    for name, (d, log_path) in LOGS.items():
        analyze(name, d, log_path)
