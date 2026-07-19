#!/usr/bin/env python3
"""Exact c_3(G) via the Laplace expansion formula.

c_3 = < G_3 + G_1 G_2 + G_1^3/6 + G_2 H_1 + G_1^2 H_1/2 + G_1 H_2 + H_3 >_mu

where:
  G_1 = d P / 12 - Q^2 / 2
  G_2 = -d^2 R / 360 + d Q P / 12 - Q^3 / 3
  G_3 = d^3 S / 20160 - d^2 Q R / 360 - d^2 P^2 / 288 + d Q^2 P / 12 - Q^4 / 4
  H_1 = -d Q / 12
  H_2 = d^2 (Q^2 / 288 - P / 1440)
  H_3 = d^3 (-Q^3/10368 + Q P / 17280 - R / 60480)

Expectations <Q^a P^b R^c S^d>_mu computed via:
  (i)   rank-2 isotropy: Q = c' |u|^2
  (ii)  rank-4 isotropy: P = (5/(d+2)) Q^2  (pointwise polynomial identity)
  (iii) rank-6 structure: R = rho_1 |u|^6 + rho_2 I_6 where rho_1 = B_3 * 15/(r(r+2)(r+4))
        For G ∈ {G_2, F_4}: <I_6>_mu = 0 rigorously, so <R Q^k>_mu = rho_1 * <|u|^6 Q^k>_mu
        (needs verification that <Q^k I_6>_mu = 0 too -- via Weyl + spherical-harmonic argument)
        For E_8: rho_2 = 0 (rank-6 1-dim)
        For E_6, E_7: <I_6>_mu != 0 -- use MC or defer
  (iv)  rank-8 structure: S = eta_1 |u|^8 + eta_2 I_8 where eta_1 similarly

For a first pass: assume <(Q^k I_{2m})>_mu = 0 for primary I_{2m} whenever applicable.
This may not be rigorous but provides a concrete rational c_3 estimate.
"""
from fractions import Fraction

ROOT_DATA = {
    'G2': {
        'r': 2, 'R_plus': 6,
        'A': 3*1 + 3*3, 'B_2': 3*1 + 3*9, 'B_3': 3*1 + 3*27, 'B_4': 3*1 + 3*81,
        'd': 14,
    },
    'F4': {
        'r': 4, 'R_plus': 24,
        'A': 12*1 + 12*2, 'B_2': 12*1 + 12*4, 'B_3': 12*1 + 12*8, 'B_4': 12*1 + 12*16,
        'd': 52,
    },
    'E6': {
        'r': 6, 'R_plus': 36,
        'A': 36*2, 'B_2': 36*4, 'B_3': 36*8, 'B_4': 36*16,
        'd': 78,
    },
    'E7': {
        'r': 7, 'R_plus': 63,
        'A': 63*2, 'B_2': 63*4, 'B_3': 63*8, 'B_4': 63*16,
        'd': 133,
    },
    'E8': {
        'r': 8, 'R_plus': 120,
        'A': 120*2, 'B_2': 120*4, 'B_3': 120*8, 'B_4': 120*16,
        'd': 248,
    },
}

def derived(G):
    rd = ROOT_DATA[G]
    d = rd['d']; r = rd['r']; Rp = rd['R_plus']
    A = rd['A']
    c_prime = Fraction(A, r)
    a = Fraction(d, 2)
    return d, r, Rp, A, c_prime, a

def moment_QK(G, K):
    """<Q^K>_mu = a(a+1)...(a+K-1) (Gamma identity, exact)."""
    _, _, _, _, _, a = derived(G)
    res = Fraction(1)
    for i in range(K):
        res *= (a + i)
    return res

def moment_u_2K(G, K):
    """<|u|^{2K}>_mu = <Q^K>/c'^K."""
    _, _, _, _, c_prime, _ = derived(G)
    return moment_QK(G, K) / (c_prime ** K)

def rho_1(G):
    """rho_1 such that sum(alpha.u)^6 = rho_1 |u|^6 + (primary I_6 contribution)."""
    _, r, _, _, _, _ = derived(G)
    rd = ROOT_DATA[G]
    B_3 = rd['B_3']
    return Fraction(B_3 * 15, r * (r+2) * (r+4))

def eta_1(G):
    """eta_1 such that sum(alpha.u)^8 = eta_1 |u|^8 + (primary I_8 contribution)."""
    _, r, _, _, _, _ = derived(G)
    rd = ROOT_DATA[G]
    B_4 = rd['B_4']
    # For |u|^8 angular average: E_Omega[hat w_1^8] = 105 / (r(r+2)(r+4)(r+6))
    # by Gaussian-sphere identity.  So eta_1 = B_4 * 105/(r(r+2)(r+4)(r+6)).
    return Fraction(B_4 * 105, r * (r+2) * (r+4) * (r+6))

def moment_P_Qk(G, k):
    """<P Q^k>_mu = (5/(d+2)) <Q^(k+2)>_mu (rank-4 isotropy, exact)."""
    d, _, _, _, _, _ = derived(G)
    return Fraction(5, d + 2) * moment_QK(G, k + 2)

def moment_Pj_Qk(G, j, k):
    """<P^j Q^k>_mu = (5/(d+2))^j <Q^(2j+k)>_mu."""
    d, _, _, _, _, _ = derived(G)
    return Fraction(5, d + 2) ** j * moment_QK(G, 2*j + k)

def moment_R_Qk(G, k, assume_I6_zero=True):
    """<R Q^k>_mu.  If we can assume <I_6 Q^k>_mu = 0:
       <R Q^k> = rho_1 <|u|^6 Q^k> = rho_1 <Q^(k+3)>/c'^3
    """
    _, _, _, _, c_prime, _ = derived(G)
    contrib_radial = rho_1(G) * moment_QK(G, k + 3) / (c_prime ** 3)
    if assume_I6_zero:
        return contrib_radial
    else:
        return None  # requires primary I_6 integration

def moment_S_Qk(G, k, assume_I8_zero=True):
    """<S Q^k>_mu with assumption <I_8 Q^k>_mu = 0."""
    _, _, _, _, c_prime, _ = derived(G)
    contrib_radial = eta_1(G) * moment_QK(G, k + 4) / (c_prime ** 4)
    if assume_I8_zero:
        return contrib_radial
    else:
        return None

def compute_c3(G, assume_zeros=True):
    """Compute c_3(G) as an exact rational, under simplifying assumption that
    primary invariant expectations <I_6 Q^k>, <I_8 Q^k> vanish.

    Notes on rigor:
    - For G_2, F_4: <I_6>_mu = 0 rigorously (via Wick / polynomial identity).
    - For E_8: rho_2 = 0 (rank-6 1-dim), so only primary I_8 enters, not I_6.
    - For E_6, E_7: <I_6>_mu != 0, so this assumption is approximate.
    - <I_6 Q^k>_mu for k > 0: may or may not vanish; not rigorously known.
    """
    d, _, _, _, _, _ = derived(G)
    d_f = Fraction(d)

    # === G_3 terms ===
    # d^3 S/20160
    t_S = d_f**3 / 20160 * moment_S_Qk(G, 0, assume_I8_zero=assume_zeros)
    # -d^2 QR/360
    t_QR = -d_f**2 / 360 * moment_R_Qk(G, 1, assume_I6_zero=assume_zeros)
    # -d^2 P^2/288
    t_P2 = -d_f**2 / 288 * moment_Pj_Qk(G, 2, 0)
    # d Q^2 P/12
    t_Q2P = d_f / 12 * moment_P_Qk(G, 2)
    # -Q^4/4
    t_Q4 = -Fraction(1, 4) * moment_QK(G, 4)
    G_3_avg = t_S + t_QR + t_P2 + t_Q2P + t_Q4

    # === G_1 G_2 expanded ===
    # G_1 G_2 = (dP/12 - Q^2/2)(-d^2 R/360 + dQP/12 - Q^3/3)
    # Term by term:
    # (dP/12)(-d^2R/360) = -d^3 PR/4320
    # (dP/12)(dQP/12) = d^2 QP^2/144
    # (dP/12)(-Q^3/3) = -dPQ^3/36
    # (-Q^2/2)(-d^2 R/360) = d^2 Q^2 R/720
    # (-Q^2/2)(dQP/12) = -dQ^3 P/24
    # (-Q^2/2)(-Q^3/3) = Q^5/6
    t = Fraction(0)
    t += -d_f**3 / 4320 * Fraction(5, d+2) * moment_R_Qk(G, 1, assume_I6_zero=assume_zeros)  # PR = 5/(d+2) Q^2 R? No wait PR is product, not rank-4 reducible
    # Actually PR = P * R, and we use P = (5/(d+2)) Q^2 pointwise. So <PR> = (5/(d+2)) <Q^2 R>.
    t += d_f**2 / 144 * Fraction(5, d+2)**2 * moment_QK(G, 5)  # QP^2 = (5/(d+2))^2 Q^5 via P = ... Q^2
    # Wait, P^2 = (5/(d+2))^2 Q^4 pointwise. So QP^2 = Q * (5/(d+2))^2 Q^4 = (5/(d+2))^2 Q^5.
    t += -d_f / 36 * moment_P_Qk(G, 3)  # PQ^3 = (5/(d+2)) Q^5 via P = ... Q^2. so PQ^3 = (5/(d+2)) Q * Q^3 = no wait PQ^3 = P * Q^3 = (5/(d+2)) Q^2 * Q^3 = (5/(d+2)) Q^5
    # Hmm, moment_P_Qk(G, 3) = (5/(d+2)) <Q^5>. yes.
    t += d_f**2 / 720 * moment_R_Qk(G, 2, assume_I6_zero=assume_zeros)  # Q^2 R
    t += -d_f / 24 * moment_P_Qk(G, 3)  # Q^3 P = same as above
    t += Fraction(1, 6) * moment_QK(G, 5)  # Q^5
    # Consolidate the t_QP^2 computation: (dP/12)(dQP/12) = d^2 Q P^2/144. <QP^2> = (5/(d+2))^2 <Q^5> ✓
    # First term: (dP/12)(-d^2 R/360) = -d^3 P R / 4320. <PR> via P = (5/(d+2)) Q^2 : <PR> = (5/(d+2)) <Q^2 R>.
    # Rewrite:
    t_G1G2 = (-d_f**3 / 4320 * Fraction(5, d+2) * moment_R_Qk(G, 2, assume_I6_zero=assume_zeros)
              + d_f**2 / 144 * Fraction(5, d+2)**2 * moment_QK(G, 5)
              - d_f / 36 * Fraction(5, d+2) * moment_QK(G, 5)
              + d_f**2 / 720 * moment_R_Qk(G, 2, assume_I6_zero=assume_zeros)
              - d_f / 24 * Fraction(5, d+2) * moment_QK(G, 5)
              + Fraction(1, 6) * moment_QK(G, 5))

    # === G_1^3/6 ===
    # (dP/12 - Q^2/2)^3/6 = (A-B)^3/6 = (A^3 - 3A^2 B + 3AB^2 - B^3)/6
    # with A = dP/12, B = Q^2/2
    # A^3 = d^3 P^3/1728 = d^3 (5/(d+2))^3 Q^6 / 1728
    # A^2 B = d^2 P^2 Q^2 / 288 = d^2 (5/(d+2))^2 Q^6 / 288
    # A B^2 = dP Q^4 / 48 = d (5/(d+2)) Q^6 / 48
    # B^3 = Q^6 / 8
    t_G1_cube = (d_f**3 * Fraction(5, d+2)**3 / 1728 * moment_QK(G, 6)
                 - 3 * d_f**2 * Fraction(5, d+2)**2 / 288 * moment_QK(G, 6)
                 + 3 * d_f * Fraction(5, d+2) / 48 * moment_QK(G, 6)
                 - Fraction(1, 8) * moment_QK(G, 6)) / 6

    # === G_2 H_1 ===
    # G_2 H_1 = (-d^2 R/360 + dQP/12 - Q^3/3)(-dQ/12)
    # = d^3 QR/(360*12) - d^2 Q^2 P/(12*12) + d Q^4/(3*12)
    t_G2H1 = (d_f**3 / (360 * 12) * moment_R_Qk(G, 1, assume_I6_zero=assume_zeros)
              - d_f**2 / 144 * moment_P_Qk(G, 2)
              + d_f / 36 * moment_QK(G, 4))

    # === G_1^2 H_1 / 2 ===
    # G_1^2 = (dP/12 - Q^2/2)^2 = d^2 P^2/144 - dPQ^2/12 + Q^4/4
    # G_1^2 H_1 = (above) * (-dQ/12)
    # = -d^3 P^2 Q/(144*12) + d^2 P Q^3/(12*12) - d Q^5/(4*12)
    # Divide by 2:
    t_G1sq_H1 = (-d_f**3 / (144*12) * Fraction(5, d+2)**2 * moment_QK(G, 5)
                 + d_f**2 / 144 * Fraction(5, d+2) * moment_QK(G, 5)
                 - d_f / 48 * moment_QK(G, 5)) / 2

    # === G_1 H_2 ===
    # G_1 H_2 = (dP/12 - Q^2/2) * d^2 (Q^2/288 - P/1440)
    # = d^3 PQ^2/(12*288) - d^3 P^2/(12*1440) - d^2 Q^4/(2*288) + d^2 PQ^2/(2*1440)
    t_G1H2 = (d_f**3 / (12*288) * moment_P_Qk(G, 2)
              - d_f**3 / (12*1440) * Fraction(5, d+2)**2 * moment_QK(G, 4)
              - d_f**2 / (2*288) * moment_QK(G, 4)
              + d_f**2 / (2*1440) * moment_P_Qk(G, 2))

    # === H_3 ===
    # H_3 = d^3(-Q^3/10368 + QP/17280 - R/60480)
    t_H3 = (-d_f**3 / 10368 * moment_QK(G, 3)
            + d_f**3 / 17280 * moment_P_Qk(G, 1)
            - d_f**3 / 60480 * moment_R_Qk(G, 0, assume_I6_zero=assume_zeros))

    c_3 = G_3_avg + t_G1G2 + t_G1_cube + t_G2H1 + t_G1sq_H1 + t_G1H2 + t_H3
    return c_3

def main():
    print(f"{'G':<4}{'d':>5}{'c_3 (assuming primary invariants vanish)':>50}")
    print('-' * 65)
    for G in ['G2', 'F4', 'E6', 'E7', 'E8']:
        c_3 = compute_c3(G, assume_zeros=True)
        c_3_f = float(c_3)
        print(f"{G:<4}{ROOT_DATA[G]['d']:>5}{str(c_3)[:30]:>30}  =  {c_3_f:.4e}")

    # Comparison to empirical
    print("\nEmpirical |R_3| * n^3 at top of character-ring range:")
    emp = {'G2': 36, 'F4': 3654, 'E6': 13670, 'E7': 83230, 'E8': 660000}
    for G in ['G2', 'F4', 'E6', 'E7', 'E8']:
        c_3 = compute_c3(G)
        print(f"  {G}: c_3_computed = {float(c_3):.3e}, |R_3| n^3 empirical = {emp[G]}, ratio = {float(c_3)/emp[G]:.3f}")

if __name__ == '__main__':
    main()
