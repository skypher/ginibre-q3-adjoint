#!/usr/bin/env python3
"""Rigorous upper bound on |c_3(G)| for each family via root-data triangle
inequality + Gamma-integral identity.

Derivation:
  m_n / tilde_m_n = < exp(sum_k G_k / n^k) * (1 + sum_k H_k / n^k) >_mu
with G_k, H_k explicit polynomials in u.  Expanding:
  c_3 = < G_3 + G_1 G_2 + G_1^3/6 + G_2 H_1 + G_1^2 H_1/2 + G_1 H_2 + H_3 >_mu

Each term is a polynomial in Q, P, R, S = Q_1, Q_2, Q_3, Q_4 where
Q_k(u) = sum_{alpha in R_+}(alpha . u)^{2k}, scaled by explicit d-powers.

Upper bound via triangle inequality on each term + pointwise
  |Q_k(u)| <= B_k * |u|^{2k}   with B_k := sum |alpha|^{2k}
and Gamma identity
  <|u|^{2K}>_mu = Gamma(a+K) / (c'^K Gamma(a))
                = a(a+1)(a+2)...(a+K-1) / c'^K

For POINTWISE-VALID bounds, we use absolute values inside every monomial.
(Note: the sign structure is ignored -- gives upper bound only.)

Similarly for H_k: expand h(t) = (exp(-Q/12 - P/1440 - R/60480 - ...))

Define the full c_3 as a sum of terms, compute triangle bound on each,
output per-family bound on |c_3(G)|.
"""
from fractions import Fraction

# Root data for each family
# A = sum |alpha|^2, B = B_2 = sum |alpha|^4, B_3 = sum |alpha|^6, B_4 = sum |alpha|^8
# for simply-laced (E_6, E_7, E_8), all |alpha|^2 = 2, so B_k = |R_+| * 2^k
# for G_2, F_4 with long/short: compute from counts
ROOT_DATA = {
    'G2': {  # rank 2, 3 short (|alpha|^2 = 1) + 3 long (|alpha|^2 = 3)
        'r': 2, 'R_plus': 6,
        'A': 3*1 + 3*3, 'B_2': 3*1 + 3*9, 'B_3': 3*1 + 3*27, 'B_4': 3*1 + 3*81,
        'd': 14,
    },
    'F4': {  # rank 4, 12 short (|alpha|^2 = 1) + 12 long (|alpha|^2 = 2)
        'r': 4, 'R_plus': 24,
        'A': 12*1 + 12*2, 'B_2': 12*1 + 12*4, 'B_3': 12*1 + 12*8, 'B_4': 12*1 + 12*16,
        'd': 52,
    },
    'E6': {  # rank 6, 36 roots, all |alpha|^2 = 2
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

def cp(G):
    """c' = A/r (Gram isotropy constant)."""
    rd = ROOT_DATA[G]
    return Fraction(rd['A'], rd['r'])

def a_param(G):
    """Gamma parameter a = d/2 = |R_+| + r/2."""
    rd = ROOT_DATA[G]
    return Fraction(rd['d'], 2)

def moment_u2K(G, K):
    """<|u|^{2K}>_mu = a(a+1)...(a+K-1) / c'^K."""
    a = a_param(G)
    c_prime = cp(G)
    num = Fraction(1)
    for i in range(K):
        num *= (a + i)
    return num / (c_prime ** K)

def bound_QaPbRcSd(G, a, b, c, d_exp):
    """Upper bound on |<Q^a P^b R^c S^{d_exp}>_mu| via triangle:
       <|Q^a P^b R^c S^{d_exp}|>_mu <= B_1^a B_2^b B_3^c B_4^{d_exp} * <|u|^{2(a+2b+3c+4d)}>_mu
    """
    rd = ROOT_DATA[G]
    B = [rd['A'], rd['B_2'], rd['B_3'], rd['B_4']]  # B_1, B_2, B_3, B_4
    K = a + 2*b + 3*c + 4*d_exp
    prefactor = Fraction(B[0])**a * Fraction(B[1])**b * Fraction(B[2])**c * Fraction(B[3])**d_exp
    return prefactor * moment_u2K(G, K)

def c3_upper_bound(G):
    """Upper bound on |c_3(G)| by triangle-inequality over all terms.

    c_3 = < G_3 + G_1 G_2 + G_1^3/6 + G_2 H_1 + G_1^2 H_1/2 + G_1 H_2 + H_3 >_mu

    where:
      G_1 = dP/12 - Q^2/2
      G_2 = -d^2 R/360 + d Q P/12 - Q^3/3
      G_3 = d^3 S/20160 - d^2 Q R/360 - d^2 P^2/288 + d Q^2 P/12 - Q^4/4
      H_1 = -d Q/12
      H_2 = d^2 Q^2/288 - d^2 P/1440
      H_3 = -d^3 Q^3/10368 + d^3 Q P/17280 - d^3 R/60480 (approximate from sin expansion)

    For rigor we use ABSOLUTE VALUE on each monomial and triangle inequality.
    """
    rd = ROOT_DATA[G]
    d = rd['d']
    d_f = Fraction(d)
    total = Fraction(0)

    # G_3 terms: d^3 S/20160, d^2 QR/360, d^2 P^2/288, d Q^2 P/12, Q^4/4
    total += abs(d_f**3) / 20160 * bound_QaPbRcSd(G, 0, 0, 0, 1)  # S
    total += abs(d_f**2) / 360   * bound_QaPbRcSd(G, 1, 0, 1, 0)  # QR
    total += abs(d_f**2) / 288   * bound_QaPbRcSd(G, 0, 2, 0, 0)  # P^2
    total += abs(d_f) / 12       * bound_QaPbRcSd(G, 2, 1, 0, 0)  # Q^2 P
    total += Fraction(1, 4)      * bound_QaPbRcSd(G, 4, 0, 0, 0)  # Q^4

    # G_1 G_2 = (dP/12 - Q^2/2)(-d^2R/360 + dQP/12 - Q^3/3)
    # Expand: dP/12 * -d^2R/360 + dP/12 * dQP/12 + dP/12 * -Q^3/3
    #        + -Q^2/2 * -d^2R/360 + -Q^2/2 * dQP/12 + -Q^2/2 * -Q^3/3
    # In absolute value, sum of |coef| * |monomial|.
    total += abs(d_f**3) / (12 * 360) * bound_QaPbRcSd(G, 0, 1, 1, 0)  # PR
    total += abs(d_f**2) / (12 * 12)  * bound_QaPbRcSd(G, 1, 2, 0, 0)  # QP^2
    total += abs(d_f) / (12 * 3)      * bound_QaPbRcSd(G, 3, 1, 0, 0)  # Q^3 P
    total += abs(d_f**2) / (2 * 360)  * bound_QaPbRcSd(G, 2, 0, 1, 0)  # Q^2 R
    total += abs(d_f) / (2 * 12)      * bound_QaPbRcSd(G, 3, 1, 0, 0)  # Q^3 P
    total += Fraction(1, 6)           * bound_QaPbRcSd(G, 5, 0, 0, 0)  # Q^5

    # G_1^3/6 = (dP/12 - Q^2/2)^3 / 6 = (d^3 P^3/1728 - 3 d^2 P^2 Q^2/(144 * 2) + 3 dP Q^4/(12 * 4) - Q^6/8)/6
    # Let me just expand: (A - B)^3 = A^3 - 3 A^2 B + 3 A B^2 - B^3
    # with A = dP/12, B = Q^2/2
    # = d^3 P^3 / 1728 - 3 d^2 P^2 Q^2 / (144*2) + 3 dP Q^4 / (12*4) - Q^6 / 8
    # Divide by 6:
    total += abs(d_f**3) / (1728 * 6)   * bound_QaPbRcSd(G, 0, 3, 0, 0)  # P^3
    total += abs(d_f**2) * 3 / (288 * 6) * bound_QaPbRcSd(G, 2, 2, 0, 0)  # Q^2 P^2
    total += abs(d_f) * 3 / (48 * 6)     * bound_QaPbRcSd(G, 4, 1, 0, 0)  # Q^4 P
    total += Fraction(1, 48)             * bound_QaPbRcSd(G, 6, 0, 0, 0)  # Q^6

    # G_2 H_1 = (-d^2 R/360 + dQP/12 - Q^3/3)(-dQ/12)
    # = d^3 Q R/(360 * 12) - d^2 Q^2 P/(12*12) + d Q^4/(3*12)
    total += abs(d_f**3) / (360 * 12)  * bound_QaPbRcSd(G, 1, 0, 1, 0)  # QR
    total += abs(d_f**2) / (12 * 12)   * bound_QaPbRcSd(G, 2, 1, 0, 0)  # Q^2 P
    total += abs(d_f) / (3 * 12)       * bound_QaPbRcSd(G, 4, 0, 0, 0)  # Q^4

    # G_1^2 H_1 / 2 = (dP/12 - Q^2/2)^2 * (-dQ/12) / 2
    # = -(dQ/12)*(d^2 P^2/144 - 2*(dP/12)(Q^2/2) + Q^4/4)/2
    # = -(dQ/24)*(d^2 P^2/144 - d Q^2 P/12 + Q^4/4)
    # = -d^3 Q P^2/(24*144) + d^2 Q^3 P/(24*12) - d Q^5/(24*4)
    total += abs(d_f**3) / (24 * 144)  * bound_QaPbRcSd(G, 1, 2, 0, 0)  # QP^2
    total += abs(d_f**2) / (24 * 12)   * bound_QaPbRcSd(G, 3, 1, 0, 0)  # Q^3 P
    total += abs(d_f) / (24 * 4)       * bound_QaPbRcSd(G, 5, 0, 0, 0)  # Q^5

    # G_1 H_2 = (dP/12 - Q^2/2)(d^2 Q^2/288 - d^2 P/1440)
    # = d^3 Q^2 P/(12 * 288) - d^3 P^2/(12 * 1440) - d^2 Q^4/(2 * 288) + d^2 Q^2 P/(2 * 1440)
    total += abs(d_f**3) / (12 * 288)  * bound_QaPbRcSd(G, 2, 1, 0, 0)  # Q^2 P
    total += abs(d_f**3) / (12 * 1440) * bound_QaPbRcSd(G, 0, 2, 0, 0)  # P^2
    total += abs(d_f**2) / (2 * 288)   * bound_QaPbRcSd(G, 4, 0, 0, 0)  # Q^4
    total += abs(d_f**2) / (2 * 1440)  * bound_QaPbRcSd(G, 2, 1, 0, 0)  # Q^2 P

    # H_3: from expansion of h(t) = exp(-Q/12 - P/1440 - R/60480 - ...)
    # h_3 coefficient in 1/n^3: from third-order Taylor of exponential.
    # log h = -Q/12 - P/1440 - R/60480 - ...
    # h = exp(log h) = 1 + log h + (log h)^2/2 + (log h)^3/6 + ...
    # To order 1/n^3 (after rescaling, Q scales as d/n, P as (d/n)^2, R as (d/n)^3):
    #   H_1 = -dQ/12
    #   H_2 = d^2 Q^2/288 - d^2 P/1440 (from (log h)^2/2 Q^2 term + linear R-like terms)
    #   H_3 = d^3(-Q^3/10368 + QP/17280 - R/60480 + ...)
    # Coefficient computation:
    #   (log h)^3/6 Q^3 term: (-1/12)^3 Q^3 / 6 = -Q^3/(12^3 * 6) = -Q^3/10368
    #   (log h)^2/2 with one Q one P: 2 * (-1/12)(-1/1440) Q P / 2 = QP/(12*1440) = QP/17280
    #   (log h) R term: -R/60480
    # After scaling:
    total += abs(d_f**3) / 10368  * bound_QaPbRcSd(G, 3, 0, 0, 0)  # Q^3
    total += abs(d_f**3) / 17280  * bound_QaPbRcSd(G, 1, 1, 0, 0)  # QP
    total += abs(d_f**3) / 60480  * bound_QaPbRcSd(G, 0, 0, 1, 0)  # R

    return total

def main():
    print(f"{'G':<4}{'d':>5}{'|c_3| <= (rigorous bound)':>35}{'(4 C_3)^(1/3)':>20}{'n_max':>10}{'close?':>8}")
    n_max_data = {'G2': 195, 'F4': 165, 'E6': 76, 'E7': 66, 'E8': 66}
    for G in ['G2', 'F4', 'E6', 'E7', 'E8']:
        bound = c3_upper_bound(G)
        bound_f = float(bound)
        threshold = (4 * bound_f) ** (1/3)
        n_max = n_max_data[G]
        ok = '✓' if threshold <= n_max else '✗'
        print(f"{G:<4}{ROOT_DATA[G]['d']:>5}{bound_f:>35.4e}{threshold:>20.2f}{n_max:>10}{ok:>8}")

if __name__ == '__main__':
    main()
