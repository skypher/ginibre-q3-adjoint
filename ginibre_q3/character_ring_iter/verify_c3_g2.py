#!/usr/bin/env python3
"""Verify the exact c_3(G_2) = -3080875/3072 (from closed-form computation)
against character-ring moments, with c_1 and c_2 FIXED at their rigorous
closed-form values.

If the closed-form c_3 is correct, then
   n^3 (m_n / tilde_m_n - 1 - c_1/n - c_2/n^2)
should converge to c_3 as n -> infinity.
"""
import re, math
import os
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))

def parse(path):
    moms = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'm_(\d+) = (\d+)', line)
            if m:
                moms[int(m.group(1))] = int(m.group(2))
    return moms

def main():
    # G_2: rigorous values
    d = 14
    c_1 = -Fraction(d * (d + 4), 16)  # -63/4
    c_2 = Fraction(165851, 1152)
    c_3_predicted = -Fraction(3080875, 3072)

    path = os.path.join(HERE, 'logs', 'g2_200.log')
    moms = parse(path)
    K = max(moms.keys())

    print(f"G_2: d = {d}, c_1 = {c_1} = {float(c_1)}")
    print(f"     c_2 = {c_2} = {float(c_2)}")
    print(f"     c_3 (predicted from closed form) = {c_3_predicted} = {float(c_3_predicted)}")
    print()
    print(f"Check: n^3 (m_n/tilde_m_n - 1 - c_1/n - c_2/n^2) -> c_3 as n -> inf")
    print()
    print(f"{'n':>5}{'n^3 * (ratio - 1 - c1/n - c2/n^2)':>42}{'diff from c_3_predicted':>30}")

    # tilde_m_n = const * d^n * n^{-d/2}
    # We don't know the const exactly, but we can eliminate it by comparing
    # ratios.  Instead, compute log(m_n / (d^(n+d/2) n^(-d/2))) - log(const)
    # and fit.

    # Simpler: compute log(m_n / (d^n n^{-d/2})) - reference_at_n_max
    # and then fit residual.  But we need absolute c_3, so we need const.
    #
    # Approach: write log(m_n/tilde_m_n) = log(const_ratio) + log(1 + c_1/n + c_2/n^2 + c_3/n^3 + ...)
    # Use LARGEST n as reference point and FIT const_ratio from c_1, c_2, c_3 assumption.

    # Actually let's just use: ratio := m_n / m_{n+2} * some scaling to cancel const.
    # Or: use CONSECUTIVE RATIOS.
    # Cleanest: compute log(m_n), subtract log(d^{n+d/2} n^{-d/2}) + log(I_0/|W|)
    # and fit log const from large-n asymptotic.

    # Let's compute log(m_n) - n log d + (d/2) log(n) and plot what it converges to.
    # That's log(const) = log(I_0/|W|) + log(1 + c_1/n + c_2/n^2 + ...)
    # For large n, this approaches log(I_0/|W|).

    const_est = None
    ns = sorted(moms.keys())
    ns = [n for n in ns if n >= 50 and n in moms and moms[n] > 0]

    # Use the LARGEST n to estimate log(I_0/|W|) from:
    # log(const) = log(m_n) - n log(d) + (d/2) log(n) - log(1 + c_1/n + c_2/n^2 + c_3/n^3)
    n_ref = max(ns)
    log_m_ref = mp.log(moms[n_ref])
    c1_mp = mp.mpf(c_1.numerator) / mp.mpf(c_1.denominator)
    c2_mp = mp.mpf(c_2.numerator) / mp.mpf(c_2.denominator)
    c3_mp = mp.mpf(c_3_predicted.numerator) / mp.mpf(c_3_predicted.denominator)
    correction_ref = mp.log(1 + c1_mp/n_ref + c2_mp/n_ref**2 + c3_mp/mp.mpf(n_ref)**3)
    log_const = log_m_ref - n_ref * mp.log(d) + (d/2) * mp.log(n_ref) - correction_ref
    print(f"  log(const) estimated from n_ref = {n_ref}: {float(log_const):.6f}")
    print()

    # Now check: for each n, compute (m_n/tilde_m_n_absolute) - 1 - c_1/n - c_2/n^2
    # where tilde_m_n_absolute = exp(log_const) * d^n * n^{-d/2}
    for n in ns[-15::2]:
        log_tilde = log_const + n * mp.log(d) - (d/2) * mp.log(n)
        ratio = mp.exp(mp.log(moms[n]) - log_tilde)
        residual = ratio - 1 - c1_mp/n - c2_mp/mp.mpf(n)**2
        n3_res = float(residual * mp.mpf(n)**3)
        diff = n3_res - float(c_3_predicted)
        print(f"{n:>5}{n3_res:>42.6f}{diff:>30.6f}")

if __name__ == '__main__':
    main()
