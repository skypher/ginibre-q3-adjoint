#!/usr/bin/env python3
"""Rigorous c_3(F_4) including primary-invariant contributions via exact Wick.

Findings:
  <I_6(F_4) Q^k>_mu = 0 exactly for k = 0, 1, 2, 3 (verified by Wick).
  <I_8(F_4) Q^k>_mu != 0 in general.

c_3 formula decomposes all the way to <Q^a P^b R^c S^d> expectations,
each of which has the form
   (coeff) * <Q^K>_mu + (primary-invariant contribution)
For R-involving terms: primary contributions = 0 for F_4.
For S-involving terms (appears only in G_3): primary contribution via
   <S_primary Q^k>_mu  computed exactly by Wick.
"""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from c2_exact import (
    f4_positive_roots, linear_form_squared_poly, linear_form_power_poly,
    poly_multiply_nd, wick_int_nd
)

def main():
    roots = f4_positive_roots()
    n_vars = 4
    d = 52
    c_prime = Fraction(9)
    a = Fraction(d, 2)

    # Set up polynomials
    prod_poly = {tuple(0 for _ in range(n_vars)): Fraction(1)}
    for alpha in roots:
        sq = linear_form_squared_poly(alpha)
        prod_poly = poly_multiply_nd(prod_poly, sq)

    # bar_S = sum(alpha.w)^8
    barS_poly = {}
    for alpha in roots:
        p8 = linear_form_power_poly(alpha, 8)
        for k, c in p8.items():
            barS_poly[k] = barS_poly.get(k, Fraction(0)) + c
    barS_poly = {k: v for k, v in barS_poly.items() if v != 0}

    # |w|^8
    w_sq = {}
    for i in range(n_vars):
        key = tuple(2 if k == i else 0 for k in range(n_vars))
        w_sq[key] = Fraction(1)
    w4 = poly_multiply_nd(w_sq, w_sq)
    w8 = poly_multiply_nd(w4, w4)

    eta_1 = Fraction(357, 32)

    # S_primary(w) = bar_S(w) - eta_1 |w|^8
    S_primary = dict(barS_poly)
    for k, c in w8.items():
        new = S_primary.get(k, Fraction(0)) - eta_1 * c
        if new != 0: S_primary[k] = new
        elif k in S_primary: del S_primary[k]

    # Q(w) = c' |w|^2 (but in w-coords where <|w|^{2k}>_mu = a(a+1)...)
    # Actually in W-coords, Q(w) = sum(alpha.w)^2 which is c'|w|^2 for rank-2 isotropy.
    # For F_4: <Q^k>_mu = <(c'|w|^2)^k>_mu = c'^k * a(a+1)...(a+k-1)
    # But we want <Q^k> in paper's u-coords = a(a+1)...(a+k-1).
    # Since u = w/sqrt(c'), Q(u) = |w|^2 in w-coords (since Q(u) = c'|u|^2 and |u|^2 = |w|^2/c').
    # So in w-coords, <Q^k>_mu = <|w|^{2k}>_mu^(w) = a(a+1)...(a+k-1).  Consistent with u-coords!

    # Proceeding: compute <S_primary(u)>, <S_primary Q>, <S_primary Q^2>, <S_primary Q^3>
    # In w-coords: <S_primary(w) * Q(u)^k>_mu^(w) = <S_primary(w) * (|w|^2)^k>_mu^(w)
    # Hmm, but S_primary(u) vs S_primary(w) differ by scale.
    # Let me just compute <S_primary(w) * |w|^{2k}>_mu^(w) which is what I need to plug in.
    #
    # The c_3 contribution from <S>_mu in paper's u-coords is d^3 <S(u)>_mu / 20160.
    # <S(u)>_mu = <eta_1 |u|^8 + S_primary(u)>_mu = eta_1 a(a+1)(a+2)(a+3)/c'^4 + <S_primary(u)>
    # And <S_primary(u)>_mu = <S_primary(w)>_mu / c'^4.  (Since S_primary degree 8 and w = u sqrt(c').)
    # So contribution: d^3/20160 * (eta_1 <|u|^8> + <S_primary(w)>/c'^4)
    #
    # For a_param style, in u-coords, <|u|^{2k}>_mu = a(a+1)...(a+k-1)/c'^k.
    # So in u-coords: <S>_mu is already in "u-units".

    # Simplest: compute <S(u)>_mu directly.
    # In w-coords: <S(w)>_mu^(w) = direct Wick value = 6265350
    # S(u) = S(w)/c'^4 (since (alpha.u)^8 = (alpha.w)^8/c'^4)
    # So <S(u)>_mu = <S(w)>_mu^(w) / c'^4 = 6265350/6561
    #
    # Hmm wait. u = w/sqrt(c'), so (alpha.u) = (alpha.w)/sqrt(c'), (alpha.u)^8 = (alpha.w)^8/c'^4.
    # S(u) = sum(alpha.u)^8 = S(w)/c'^4.
    # <S(u)>_mu: under measure mu (same physical measure, just different coords), expectation is
    # <S(w)>_mu^(w) / c'^4. ✓

    sigma_sq = Fraction(1, 2)
    Z = wick_int_nd(prod_poly, sigma_sq)

    # Compute in w-coords all needed <S Q^k>
    # Q(w) = c' |w|^2 in rank-2 isotropic case, so for our product:
    # <S(w) Q(w)^k>_mu = c'^k <S(w) |w|^{2k}>_mu
    # and similarly <|w|^{2j}>_mu = a(a+1)...(a+j-1)

    # Compute <S Q^k>_mu in u-coords (which is what c_3 formula uses):
    # <S Q^k>_mu_u = <S(u) Q(u)^k>_mu_u
    # Translation to w-coords: S(u) Q(u)^k = (S(w)/c'^4)(|w|^2)^k = S(w) |w|^{2k}/c'^4
    # (note: Q(u) = |w|^2 in w-coords, so Q(u)^k = |w|^{2k})
    # So <S(u) Q(u)^k>_mu = <S(w) |w|^{2k}>_mu^(w) / c'^4

    print(f"F_4: d = {d}, c' = {c_prime}, a = {a}, eta_1 = {eta_1}")

    # Compute <bar_S(w) |w|^{2k}>_mu^(w) via Wick for k = 0, 1, 2
    for k in range(3):
        # Build |w|^{2k}
        wk_poly = {tuple(0 for _ in range(n_vars)): Fraction(1)}
        for _ in range(k):
            wk_poly = poly_multiply_nd(wk_poly, w_sq)
        # <bar_S * |w|^{2k}>_mu
        num_full = wick_int_nd(poly_multiply_nd(poly_multiply_nd(prod_poly, barS_poly), wk_poly), sigma_sq)
        full = Fraction(num_full, Z)
        # Convert to <S(u) Q(u)^k> in u-coords:
        # <S(w) |w|^{2k}>_mu / c'^4 = <S(u)|w|^{2k}_w>/c'^4
        # But |w|^{2k} = (c'|u|^2)^k = c'^k |u|^{2k} -- wait, depends on which var |w| is.
        # Ugh. Let me just keep it w-native and convert at the end.
        # <S(u) Q(u)^k>_u-coords = <S(w)/c'^4 * |w|^{2k}>_w = <S(w) |w|^{2k}>_w / c'^4
        full_u = full / c_prime**4
        # And radial-only: eta_1 <|u|^8 Q(u)^k>_u = eta_1 <|w|^{2k+8}/c'^{4}>_w * ? wait
        # <|u|^8>_u = <|u|^8>_u = <|w|^8/c'^4>_w. So = <|w|^8>_w/c'^4.
        # <|u|^8 Q(u)^k>_u = <|w|^8 |w|^{2k}/c'^{4+k}>_w... no, Q(u) = |w|^2 in w-coords
        # Oh wait Q(u) = c'|u|^2. In w-coords, Q(u) = c'|w|^2/c' = |w|^2. So Q(u) = |w|^2 exactly.
        # So <|u|^8 Q(u)^k>_u = <(|w|^8/c'^4) |w|^{2k}>_w = <|w|^{2k+8}>_w / c'^4

        # Thus <S(u) Q(u)^k>_u = <S(w) |w|^{2k}>_w / c'^4
        # eta_1 <|u|^8 Q(u)^k>_u = eta_1 <|w|^{2k+8}>_w / c'^4 = eta_1 a(a+1)...(a+k+3) / c'^4 (Gamma in w)
        # <S_primary(u) Q^k> = full_u - (radial part)
        num_w_2k_plus_8 = Fraction(1)
        for i in range(k + 4):  # a(a+1)...(a+k+3)
            num_w_2k_plus_8 *= (a + i)
        radial_u = eta_1 * num_w_2k_plus_8 / c_prime**4
        sprim_u = full_u - radial_u

        print(f"  k = {k}: <S(u) Q^k>_u = {full_u}  = {float(full_u):.4f}")
        print(f"         eta_1 <|u|^8 Q^k>_u (radial) = {radial_u} = {float(radial_u):.4f}")
        print(f"         <S_primary(u) Q^k>_u (correction) = {sprim_u} = {float(sprim_u):.4f}")
    print()

    # In c_3 formula, S appears ONLY in G_3 term (coeff d^3/20160).
    # So total S-correction to c_3: d^3/20160 * <S_primary(u)>_mu
    # (<S_primary Q^k>_mu for k > 0 don't directly appear in c_3 — only <S> at k=0.)

    # Wait, does S appear in any other c_3 term?  Let me check:
    # c_3 = G_3 + G_1 G_2 + G_1^3/6 + G_2 H_1 + G_1^2 H_1/2 + G_1 H_2 + H_3
    # G_3 has S.  G_1, G_2, H_1, H_2, H_3 do not contain S (their highest Q_k is R = Q_3).
    # So S contributes only via G_3 term d^3 S/20160.

    # Compute rigorous S correction to c_3
    k = 0
    wk_poly = {tuple(0 for _ in range(n_vars)): Fraction(1)}
    num_full = wick_int_nd(poly_multiply_nd(poly_multiply_nd(prod_poly, barS_poly), wk_poly), sigma_sq)
    full = Fraction(num_full, Z)
    full_u = full / c_prime**4
    num_w_8 = a * (a+1) * (a+2) * (a+3)
    radial_u = eta_1 * num_w_8 / c_prime**4
    sprim_u_k0 = full_u - radial_u
    S_correction_to_c3 = Fraction(d**3, 20160) * sprim_u_k0
    print(f"RIGOROUS S-primary correction to c_3(F_4): {S_correction_to_c3} ≈ {float(S_correction_to_c3):.4f}")

    # Load radial-only c_3 from other script
    from c3_exact import compute_c3
    c3_radial = compute_c3('F4', assume_zeros=True)
    c3_rigorous = c3_radial + S_correction_to_c3
    print(f"\nc_3^radial(F_4) = {c3_radial} ≈ {float(c3_radial):.4f}")
    print(f"S correction = {S_correction_to_c3} ≈ {float(S_correction_to_c3):.4f}")
    print(f"c_3^rigorous(F_4) = {c3_rigorous} ≈ {float(c3_rigorous):.4f}")

    # Compare to empirical fit
    print(f"\nEmpirical fit (from verify_c3_g2.py style, using c_1, c_2 fixed): c_3(F_4) ≈ -1140851")

if __name__ == '__main__':
    main()
