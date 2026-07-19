#!/usr/bin/env python3
"""High-precision Monte Carlo / importance-sampled computation of c_2(G)
for G = E_6 and G = E_7 using the Laplace-expansion formula.

The measure is mu(dw) ~ e^{-|w|^2} prod(alpha.w)^2 dw on R^r (r = rank).
For simply-laced G (all |alpha|^2 = 2), we generate roots in 8D ambient
coords via Weyl-orbit BFS from Bourbaki simple roots, and integrate
using either:
  (a) Direct Monte Carlo with importance sampling from e^{-|w|^2}, or
  (b) Weyl-orbit reduction: <sum(alpha.w)^6>_mu = |R_+| * <(alpha_0.w)^6>_mu,
      which is a 1D marginal moment computable from the distribution of
      y := alpha_0.w under mu.

We do (a) with large sample count N and report statistical error via CLT.

Formula:
  c_2 = coef_barR * <sum(alpha.w)^6>_mu
      + coef_w6 * <|w|^6>_mu  [= Gamma: a(a+1)(a+2)]
      + coef_w4 * <|w|^4>_mu  [= Gamma: a(a+1)]
      + coef_w8 * <|w|^8>_mu  [= Gamma: a(a+1)(a+2)(a+3)]

with (from c2_exact.py derivation):
  coef_barR = -d^2 / (360 c'^3)
  coef_w6 = 5d/(12(d+2)) - 1/3 + (d+12)d/(144(d+2))
  coef_w4 = d^2 (d+1) / (288(d+2))
  coef_w8 = (d+12)^2 / (288(d+2)^2)

Only <sum(alpha.w)^6>_mu is estimated numerically; the rest are closed
form.

For simply-laced, <sum(alpha.w)^6>_mu = |R_+| <(alpha_0.w)^6>_mu.
"""
import os, sys, math
from fractions import Fraction
import numpy as np

# Optional: mpmath for high-precision 1D quadrature alternative
try:
    import mpmath as mp
    mp.mp.dps = 40
    HAVE_MPMATH = True
except ImportError:
    HAVE_MPMATH = False

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

def g2_simples():
    """G_2 simple roots in 2D basis.  Short: alpha_1 = e_1.  Long: alpha_2."""
    return np.array([
        [1.0, 0.0],           # alpha_1 short
        [-1.5, math.sqrt(3)/2],  # alpha_2 long (Bourbaki convention)
    ], dtype=np.float64)

def f4_simples():
    """F_4 simple roots in 4D."""
    return np.array([
        [0, 1, -1, 0],  # alpha_1 long (squared length 2)
        [0, 0, 1, -1],  # alpha_2 long
        [0, 0, 0, 1],   # alpha_3 short (squared length 1)
        [0.5, -0.5, -0.5, -0.5],  # alpha_4 short
    ], dtype=np.float64)

def e6_simples():
    return np.array([
        [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
    ], dtype=np.float64)

def e7_simples():
    return np.array([
        [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, 0, 0],
    ], dtype=np.float64)

def orbit_roots(simples, tol=1e-9):
    """Generate full root system via Weyl orbit BFS (float version).
    Uses quantized keys (round to tol^{-1}) to handle numerical precision."""
    def key(v):
        return tuple(round(float(x) / tol) for x in v)
    roots = {}
    for s in simples:
        roots[key(s)] = s.tolist()
    frontier = list(roots.values())
    while frontier:
        new_frontier = []
        for beta in frontier:
            b = np.array(beta)
            for alpha in simples:
                num = 2 * np.dot(b, alpha) / np.dot(alpha, alpha)
                new_vec = b - num * alpha
                k = key(new_vec)
                if k not in roots:
                    roots[k] = new_vec.tolist()
                    new_frontier.append(new_vec.tolist())
        frontier = new_frontier
    return np.array(list(roots.values()))

def positive_roots(roots):
    """Extract positive roots (first nonzero coord positive)."""
    out = []
    tol = 1e-10
    for r in roots:
        for x in r:
            if abs(x) > tol:
                if x > 0:
                    out.append(r)
                break
    return np.array(out)

def c2_formula_coeffs(d, r, A):
    """Return (coef_barR, coef_w6, coef_w4, coef_w8) as exact Fractions."""
    cp = Fraction(A, r)
    coef_barR = -Fraction(d*d, 360 * cp**3)
    coef_w6 = Fraction(5*d, 12 * (d+2)) - Fraction(1, 3) + Fraction((d+12)*d, 144*(d+2))
    coef_w4 = Fraction(d*d * (d+1), 288 * (d+2))
    coef_w8 = Fraction((d+12)**2, 288 * (d+2)**2)
    return coef_barR, coef_w6, coef_w4, coef_w8

def gamma_moment(a, k):
    res = Fraction(1)
    for i in range(k):
        res *= (a + i)
    return res

def importance_sample_barR(root_vectors, N, seed=0, batch=100000):
    """Monte Carlo estimate of <sum_alpha(alpha.w)^6>_mu with correct delta-method
    variance including the covariance between num and denom.
    """
    rng = np.random.default_rng(seed)
    R = root_vectors
    total_num = 0.0
    total_denom = 0.0
    M2_num = 0.0
    M2_denom = 0.0
    M_cross = 0.0
    count = 0
    dim = R.shape[1]
    # Also accumulate |w|^4 for sanity check against Gamma identity
    accum_w4_num = 0.0
    # Proposal: Gaussian with variance SIGMA^2 I (instead of 1/2 I).
    # Use log-weights to avoid float overflow/underflow.  For each
    # sample, log_weight = sum(log((alpha.w)^2)) - |w|^2 / 2
    # Then ratio = sum(f * exp(log_w - log_w_max)) / sum(exp(log_w - log_w_max))
    # (log-sum-exp stabilization)
    SIGMA = 1.0
    corr_factor = 1.0 - 1.0 / (2 * SIGMA * SIGMA)  # 1/2 for SIGMA=1
    log_w_all = []
    obs_all = []
    w_sq_all = []
    for _ in range(N // batch):
        W = rng.standard_normal(size=(batch, dim)) * SIGMA
        inner = W @ R.T
        w_sq = (W * W).sum(axis=1)
        # log |inner|^2 with safeguard against zero
        inner_sq = inner * inner
        # Avoid log(0): add tiny epsilon but track its impact
        mask = inner_sq > 0
        log_prod = np.where(mask, np.log(np.where(mask, inner_sq, 1.0)), -np.inf).sum(axis=1)
        log_w = log_prod - w_sq * corr_factor
        log_w_all.append(log_w)
        obs_all.append(np.sum(inner**6, axis=1))
        w_sq_all.append(w_sq)
        count += batch
    log_w_all = np.concatenate(log_w_all)
    obs_all = np.concatenate(obs_all)
    w_sq_all = np.concatenate(w_sq_all)
    # Log-sum-exp stabilize
    log_w_max = np.max(log_w_all)
    stable_w = np.exp(log_w_all - log_w_max)
    weights = stable_w  # Proportional to true weights (normalizer is exp(log_w_max))
    # Now recompute totals with stable weights
    num_vals = obs_all * weights
    denom_vals = weights
    total_num = num_vals.sum()
    total_denom = denom_vals.sum()
    M2_num = (num_vals ** 2).sum()
    M2_denom = (denom_vals ** 2).sum()
    M_cross = (num_vals * denom_vals).sum()
    accum_w4_num = (w_sq_all * w_sq_all * weights).sum()
    # Sanity check: <|w|^4>_mu should equal a(a+1)
    w4_mu = accum_w4_num / total_denom
    print(f"    [MC diag] <|w|^4>_mu = {w4_mu:.4f}  (Gamma expected = a(a+1))")
    mean_num = total_num / count
    mean_denom = total_denom / count
    var_num = M2_num / count - mean_num**2
    var_denom = M2_denom / count - mean_denom**2
    cov_nd = M_cross / count - mean_num * mean_denom
    ratio = mean_num / mean_denom
    # Delta method: Var(N/D) / count = (N/D)^2 [Var(N)/N^2 + Var(D)/D^2 - 2 Cov(N,D)/(N D)] / count
    # But proper expression needs signs; use Taylor expansion of N/D around means.
    # Variance of sample mean ratio: sigma^2 / count where
    # sigma^2 = (1/D^2) [Var(N) - 2 (N/D) Cov(N,D) + (N/D)^2 Var(D)]
    sigma_sq = (var_num - 2 * ratio * cov_nd + ratio**2 * var_denom) / (mean_denom**2)
    std_err = math.sqrt(max(sigma_sq, 0) / count)
    return ratio, std_err

def compute_c2_numerical(G_name, N=2_000_000, seed=0):
    if G_name == 'G2':
        simples = g2_simples()
        Rp, d, r, A = 6, 14, 2, 12
        expected_paper = 143.967882  # exact rational 165851/1152
    elif G_name == 'F4':
        simples = f4_simples()
        Rp, d, r, A = 24, 52, 4, 36
        expected_paper = 17287.458847736626  # exact 8401705/486
    elif G_name == 'E6':
        simples = e6_simples()
        Rp, d, r, A = 36, 78, 6, 72
        expected_paper = 82150.0
    elif G_name == 'E7':
        simples = e7_simples()
        Rp, d, r, A = 63, 133, 7, 126
        expected_paper = 563000.0
    else:
        raise ValueError(G_name)
    roots_all = orbit_roots(simples)
    pos = positive_roots(roots_all)
    print(f"  {G_name}: positive roots found = {len(pos)} (expected {Rp})")
    assert len(pos) == Rp

    coef_barR, coef_w6, coef_w4, coef_w8 = c2_formula_coeffs(d, r, A)
    a = Fraction(d, 2)
    Ew4 = gamma_moment(a, 2)
    Ew6 = gamma_moment(a, 3)
    Ew8 = gamma_moment(a, 4)

    # Closed-form contribution
    c2_closed_part = coef_w6 * Ew6 + coef_w4 * Ew4 + coef_w8 * Ew8
    c2_closed_float = float(c2_closed_part)
    print(f"  {G_name}: closed-form (Gamma) contribution to c_2 = {c2_closed_float}")

    # MC for <bar R>_mu
    print(f"  {G_name}: running MC (N = {N}, sigma^2 = 1/2 Gaussian) ...")
    barR_mean, barR_err = importance_sample_barR(pos, N=N, seed=seed)
    print(f"  {G_name}: <sum_alpha(alpha.w)^6>_mu ≈ {barR_mean:.4g} ± {barR_err:.4g}")

    # c_2 = coef_barR * <bar R> + closed-form part
    c2_numerical = float(coef_barR) * barR_mean + c2_closed_float
    c2_err = abs(float(coef_barR)) * barR_err
    c1 = -d * (d + 4) / 16
    disc = c1**2 - 2 * c2_numerical
    disc_err = 2 * c2_err

    print(f"  {G_name}: c_2 ≈ {c2_numerical:.4f} ± {c2_err:.4f}")
    print(f"  {G_name}: paper's value = {expected_paper}")
    print(f"  {G_name}: c_1^2 = {c1**2}")
    print(f"  {G_name}: c_1^2 - 2 c_2 ≈ {disc:.4f} ± {disc_err:.4f}")
    sign = 'NEGATIVE (disc safely < 0, series >= 1/2 for all n > 0)' if disc + 5 * disc_err < 0 else ('POSITIVE (disc > 0, series threshold needed)' if disc - 5 * disc_err > 0 else 'INCONCLUSIVE (margin within error)')
    print(f"  {G_name}: discriminant sign at 5-sigma: {sign}")
    return c2_numerical, c2_err, disc, disc_err

def compute_c2_multiseed(G_name, N, n_seeds=5):
    """Run MC with multiple seeds to get robust mean and std."""
    c2s = []
    for seed in range(n_seeds):
        c2, c2_err, disc, disc_err = compute_c2_numerical(G_name, N=N, seed=seed)
        c2s.append(c2)
    import numpy as np
    mean_c2 = np.mean(c2s)
    std_c2 = np.std(c2s, ddof=1)
    std_err = std_c2 / np.sqrt(n_seeds)
    return mean_c2, std_err, c2s

if __name__ == '__main__':
    N_mc = int(os.environ.get('N', '2000000'))
    groups = os.environ.get('GROUPS', 'G2,F4,E6,E7').split(',')
    multiseed = os.environ.get('MULTISEED', '0') == '1'
    n_seeds = int(os.environ.get('NSEEDS', '5'))
    if multiseed:
        print(f"=== Multi-seed MC c_2 (N = {N_mc} per seed, {n_seeds} seeds) ===")
        for G in groups:
            print(f"\n{G}:")
            mean_c2, std_err, c2s = compute_c2_multiseed(G, N=N_mc, n_seeds=n_seeds)
            print(f"  {G}: c_2 mean = {mean_c2:.4f} (stderr across seeds = {std_err:.4f})")
            print(f"  {G}: individual c_2 values: {[f'{c:.2f}' for c in c2s]}")
    else:
        print(f"=== Monte Carlo c_2 (N = {N_mc}) ===")
        for G in groups:
            print(f"\n{G}:")
            c2, c2_err, disc, disc_err = compute_c2_numerical(G, N=N_mc, seed=42)
