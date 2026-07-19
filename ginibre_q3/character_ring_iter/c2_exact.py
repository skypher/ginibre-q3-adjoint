#!/usr/bin/env python3
"""Exact rational computation of c_2(G) from Wong/Laplace expansion formula.

The second Laplace coefficient c_2(G) in m_n = tilde_m_n * (1 + c_1/n + c_2/n^2 + ...)
has the decomposition:

  c_2(G) = A_6 <|w|^6> + B_6 <I_6(w)>
        + (d^2 (d+1)/(288 (d+2))) <|w|^4>
        + ((d+12) d / (144(d+2))) <|w|^6>
        + ((d+12)^2 / (288(d+2)^2)) <|w|^8>

where:
  <|w|^(2k)> = a * (a+1) * ... * (a+k-1) with a = d/2         (Gamma identity)
  A_6 = -d^2 rho_1 / (360 c'^3) + 5 d / (12 (d+2)) - 1/3
  B_6 = -d^2 rho_2 / (360 c'^3)
  c' = A / r  (rank-2 Gram isotropy)
  rho_1, rho_2: decomposition R(u) = rho_1 (Q/c')^3 + rho_2 I_6(u)
                 (rank-6 invariant space)
  <I_6(w)>: Wick expectation under e^{-|w|^2} prod(alpha·w)^2 dw

For E_8: rank-6 space is 1-dim (no primary deg-6 invariant), so rho_2 = 0
and c_2(E_8) reduces to a pure Gamma-factor closed form.

For G_2, F_4, E_6, E_7: rank-6 space is 2-dim; B_6 != 0 and <I_6> must be
computed via Wick contraction against the weighted Gaussian measure.

This script:
  1. Verifies c_2(E_8) = 346377539/45 using only Gamma identities.
  2. Computes c_2(G_2) exactly via Wick (rank 2, 6 positive roots).
"""
from fractions import Fraction
from math import factorial

# ---------- Root data (rank, |R_+|, A = sum |alpha|^2, B = sum |alpha|^4) ----------
ROOT_DATA = {
    'G2': {'r': 2, 'R_plus': 6,  'A': Fraction(12),   'B': Fraction(30)},
    'F4': {'r': 4, 'R_plus': 24, 'A': Fraction(36),   'B': Fraction(60)},
    'E6': {'r': 6, 'R_plus': 36, 'A': Fraction(72),   'B': Fraction(144)},
    'E7': {'r': 7, 'R_plus': 63, 'A': Fraction(126),  'B': Fraction(252)},
    'E8': {'r': 8, 'R_plus': 120,'A': Fraction(240),  'B': Fraction(480)},
}

def derived(G):
    rd = ROOT_DATA[G]
    r, Rp, A, B = rd['r'], rd['R_plus'], rd['A'], rd['B']
    d = r + 2 * Rp
    a = Fraction(d, 2)
    c_prime = A / r
    # lambda in paper: T^(4)_ijkl = lambda (delta delta + delta delta + delta delta)
    # P(u) = sum(alpha.u)^4 = 3 lambda |u|^4, so lambda = B/(r(r+2))
    lam = B / (r * (r + 2))
    return r, Rp, A, B, d, a, c_prime, lam

def gamma_moment(a, k):
    """<|w|^{2k}> = a(a+1)...(a+k-1) under e^{-|w|^2} prod(alpha w)^2 dw normalized."""
    res = Fraction(1)
    for i in range(k):
        res *= (a + i)
    return res

def c1_closed(G):
    """Verify c_1(G) = -d(d+4)/16 via the paper's formula."""
    r, Rp, A, B, d, a, cp, lam = derived(G)
    # c_1 = (d/12)(E[P] - E[Q]) - (1/2) E[Q^2]
    # E[Q] = a, E[Q^2] = a(a+1)
    # E[P] = (3 lambda / c'^2) E[Q^2] = (3 lam / cp^2) a(a+1)
    EQ = a
    EQ2 = a * (a + 1)
    EP = (3 * lam / (cp * cp)) * EQ2
    c1 = Fraction(d, 12) * (EP - EQ) - Fraction(1, 2) * EQ2
    return c1, EP, EQ2

def c2_for_E8():
    """c_2(E_8): rank-6 invariant space is 1-dim, so rho_2 = 0 and
    rho_1 = 15 (from R = sum(alpha.u)^6 = 15|u|^6).  All <|w|^{2k}>
    are Gamma factors."""
    G = 'E8'
    r, Rp, A, B, d, a, cp, lam = derived(G)
    # rank-6 isotropy for E_8: R(u) = 15 |u|^6 exactly, so rho_1 = 15, rho_2 = 0.
    rho_1 = Fraction(15)
    rho_2 = Fraction(0)
    # A_6 = -d^2 rho_1 / (360 c'^3) + 5 d / (12(d+2)) - 1/3
    A_6 = -Fraction(d*d) * rho_1 / (360 * cp**3) + Fraction(5 * d, 12 * (d+2)) - Fraction(1, 3)
    B_6 = -Fraction(d*d) * rho_2 / (360 * cp**3)
    # <|w|^(2k)>
    Ew4 = gamma_moment(a, 2)  # a(a+1)
    Ew6 = gamma_moment(a, 3)  # a(a+1)(a+2)
    Ew8 = gamma_moment(a, 4)
    # H_2 coeff * <|w|^4> + G_1 H_1 coeff * <|w|^6> + (G_1^2/2) coeff * <|w|^8>
    H2_coef = Fraction(d*d * (d+1), 288 * (d+2))
    G1H1_coef = Fraction((d+12) * d, 144 * (d+2))
    G1sq_coef = Fraction((d+12)**2, 288 * (d+2)**2)
    # <I_6> is 0 since B_6 = 0 for E_8
    c2 = A_6 * Ew6 + B_6 * Fraction(0) + H2_coef * Ew4 + G1H1_coef * Ew6 + G1sq_coef * Ew8
    return c2, A_6, B_6

# ---------- Exact Wick for rank-2 (G_2) primary degree-6 invariant ----------

def wick_int(poly_coeffs, sigma_sq=Fraction(1, 2)):
    """Compute Gaussian expectation in 2D with Cov(v_i, v_j) = sigma_sq delta_ij.

    poly_coeffs: dict {(k1, k2): coefficient} for polynomial sum c_{k1,k2} v_1^{k_1} v_2^{k_2}.
    Returns exact sum of E[v_1^{k_1} v_2^{k_2}] weighted by coefficients.
    E[v_i^(2m)] = (2m-1)!! * sigma^(2m). Zero if any odd.
    """
    total = Fraction(0)
    for (k1, k2), c in poly_coeffs.items():
        if k1 % 2 or k2 % 2:
            continue
        m1 = k1 // 2
        m2 = k2 // 2
        # (2m-1)!! = prod_{i=1..m} (2i - 1)
        dd1 = Fraction(1)
        for i in range(1, m1 + 1):
            dd1 *= (2*i - 1)
        dd2 = Fraction(1)
        for i in range(1, m2 + 1):
            dd2 *= (2*i - 1)
        total += c * dd1 * dd2 * sigma_sq**(m1 + m2)
    return total

def poly_multiply(p, q):
    """Multiply two polynomials represented as dicts {(k1, k2): coef}."""
    result = {}
    for (a1, a2), ca in p.items():
        for (b1, b2), cb in q.items():
            key = (a1 + b1, a2 + b2)
            result[key] = result.get(key, Fraction(0)) + ca * cb
    # Clean zero coefficients
    result = {k: v for k, v in result.items() if v != 0}
    return result

def poly_power(p, n):
    """Raise polynomial (dict) to power n."""
    result = {(0, 0): Fraction(1)}
    for _ in range(n):
        result = poly_multiply(result, p)
    return result

def compute_G2_c2_via_wick():
    """Compute c_2(G_2) exactly: need <I_6(w)> via Wick contraction against
    the weighted Gaussian e^{-|w|^2} prod_{alpha}(alpha.w)^2 dw.

    G_2 has rank 2, 6 positive roots (3 short of length 1, 3 long of sqrt(3))
    in a specific basis.  I_6^{(0)}(w) = w_1^6 - 15 w_1^4 w_2^2 + 15 w_1^2 w_2^4 - w_2^6
    (dihedral D_6 primary deg-6 invariant, = Re((w_1+i w_2)^6)).

    Derive rho_1, rho_2 from R(u) = sum(alpha.u)^6 expansion:
      R(u) = 3/16 * (127 u_1^6 + 615 u_1^4 u_2^2 + 225 u_1^2 u_2^4 + 153 u_2^6)
      = (105/4) |u|^6 + (-39/16) I_6^{(0)}(u)
    (derived by hand, verified at u=(1,0) and u=(0,1))
    """
    G = 'G2'
    r, Rp, A, B, d, a, cp, lam = derived(G)
    # G_2: rho_1 = 105/4, rho_2 = -39/16
    rho_1 = Fraction(105, 4)
    rho_2 = Fraction(-39, 16)
    A_6 = -Fraction(d*d) * rho_1 / (360 * cp**3) + Fraction(5 * d, 12 * (d+2)) - Fraction(1, 3)
    B_6 = -Fraction(d*d) * rho_2 / (360 * cp**3)

    # Wick: compute <I_6^{(0)}(w)> under e^{-|w|^2} prod(alpha.w)^2 dw
    # Step 1: explicit positive roots of G_2 in 2D
    # Short: alpha_1 = (1,0), alpha_1+alpha_2 = (-1/2, sqrt(3)/2), 2alpha_1+alpha_2 = (1/2, sqrt(3)/2)
    # Long:  alpha_2 = (-3/2, sqrt(3)/2), 3alpha_1+alpha_2 = (3/2, sqrt(3)/2), 3alpha_1+2alpha_2 = (0, sqrt(3))
    # Use rational squared: work with (a, b) s.t. alpha = (a, b*sqrt(3)), so
    # (alpha.w)^2 = (a w_1 + b sqrt(3) w_2)^2 = a^2 w_1^2 + 2 a b sqrt(3) w_1 w_2 + 3 b^2 w_2^2
    # Products of these with different alpha will mix rational and sqrt(3) parts.
    # Overall, the PRODUCT prod(alpha.w)^2 is a polynomial in w_1, w_2 with rational coefficients
    # (since each factor has sqrt(3) symmetry and the product is symmetric under sqrt(3) -> -sqrt(3)).
    # We compute symbolically keeping (rational + sqrt(3)*rational) form.
    #
    # Actually simpler: work in complex coords z = w_1 + i w_2. Each root alpha has a complex
    # form, and (alpha.w) = Re(alpha_C \bar z) or similar. But let's just compute in reals.

    # Represent (alpha.w)^2 as a polynomial in w_1, w_2 over Q[sqrt(3)].
    # Expand each as a*w_1 + b*w_2 where coefficients may be ±sqrt(3)/2, etc.
    # We'll use the form: alpha = (p, q sqrt(3)) where p, q are rationals.
    # Then (alpha.w)^2 = p^2 w_1^2 + 2 p q sqrt(3) w_1 w_2 + 3 q^2 w_2^2
    # In basis {1, sqrt(3)}: coefficient of w_1^2 is p^2 (rational part), of w_1 w_2 is 2pq (sqrt(3) part), of w_2^2 is 3q^2 (rational).

    # Represent each polynomial factor as a pair (R_part, S_part) where R is rational part, S is coefficient of sqrt(3).
    # A polynomial is a dict {(k1, k2): (R_frac, S_frac)}.

    # Multiplication: (R1 + sqrt(3) S1)(R2 + sqrt(3) S2) = (R1*R2 + 3*S1*S2) + sqrt(3)*(R1*S2 + S1*R2)

    def mul_alg(p, q):
        """Multiply two polynomials over Q[sqrt 3]."""
        result = {}
        for (a1, a2), (R1, S1) in p.items():
            for (b1, b2), (R2, S2) in q.items():
                key = (a1+b1, a2+b2)
                R_new = R1 * R2 + 3 * S1 * S2
                S_new = R1 * S2 + S1 * R2
                if key in result:
                    R_old, S_old = result[key]
                    result[key] = (R_old + R_new, S_old + S_new)
                else:
                    result[key] = (R_new, S_new)
        return {k: (R, S) for k, (R, S) in result.items() if R != 0 or S != 0}

    # Positive roots of G_2 in (p, q sqrt(3)) form:
    roots = [
        (Fraction(1),      Fraction(0)),         # alpha_1
        (Fraction(-1, 2),  Fraction(1, 2)),      # alpha_1 + alpha_2
        (Fraction(1, 2),   Fraction(1, 2)),      # 2 alpha_1 + alpha_2
        (Fraction(-3, 2),  Fraction(1, 2)),      # alpha_2
        (Fraction(3, 2),   Fraction(1, 2)),      # 3 alpha_1 + alpha_2
        (Fraction(0),      Fraction(1)),         # 3 alpha_1 + 2 alpha_2
    ]

    # Build (alpha.w)^2 = p^2 w_1^2 + 2pq sqrt(3) w_1 w_2 + 3 q^2 w_2^2
    # Polynomial over Q[sqrt 3]: {(2,0): (p^2, 0), (1,1): (0, 2pq), (0,2): (3 q^2, 0)}
    def root_squared(p, q):
        return {
            (2, 0): (p*p, Fraction(0)),
            (1, 1): (Fraction(0), 2*p*q),
            (0, 2): (3*q*q, Fraction(0)),
        }

    # Compute prod(alpha.w)^2 as polynomial over Q[sqrt 3]
    prod_poly = {(0, 0): (Fraction(1), Fraction(0))}
    for (p, q) in roots:
        prod_poly = mul_alg(prod_poly, root_squared(p, q))

    # The final product should have S (sqrt 3) part = 0 by W-invariance.
    # Verify:
    s_parts_all_zero = all(S == 0 for _, (_, S) in prod_poly.items())
    if not s_parts_all_zero:
        print("  WARNING: prod_poly has nonzero sqrt(3) components (W-invariance check failed).")
    # Extract rational part
    prod_rat = {k: R for k, (R, _) in prod_poly.items() if R != 0}

    # Now compute <prod_rat * I_6^{(0)}> and <prod_rat> (normalizer) via Wick
    I6 = {
        (6, 0): Fraction(1),
        (4, 2): Fraction(-15),
        (2, 4): Fraction(15),
        (0, 6): Fraction(-1),
    }
    # Polynomial multiplication over Q
    num_poly = poly_multiply(prod_rat, I6)
    denom_poly = prod_rat

    # Wick: sigma^2 = 1/2 so that Var(v_i) = 1/2 under e^{-|v|^2}
    sigma_sq = Fraction(1, 2)
    E_prod_I6 = wick_int(num_poly, sigma_sq=sigma_sq)
    E_prod = wick_int(denom_poly, sigma_sq=sigma_sq)
    # <I_6>_mu = E_prod_I6 / E_prod
    E_I6 = E_prod_I6 / E_prod

    # Also compute normalizer I_0' = E_prod * (pi^{r/2}) — but we only need ratios,
    # so E_I6 as Fraction is fine.

    # Now assemble c_2(G_2)
    Ew4 = gamma_moment(a, 2)
    Ew6 = gamma_moment(a, 3)
    Ew8 = gamma_moment(a, 4)
    H2_coef = Fraction(d*d * (d+1), 288 * (d+2))
    G1H1_coef = Fraction((d+12) * d, 144 * (d+2))
    G1sq_coef = Fraction((d+12)**2, 288 * (d+2)**2)
    c2 = A_6 * Ew6 + B_6 * E_I6 + H2_coef * Ew4 + G1H1_coef * Ew6 + G1sq_coef * Ew8

    return c2, {
        'A_6': A_6, 'B_6': B_6,
        'Ew6': Ew6, 'Ew4': Ew4, 'Ew8': Ew8,
        'E_I6': E_I6,
        'H2_coef': H2_coef, 'G1H1_coef': G1H1_coef, 'G1sq_coef': G1sq_coef,
    }

def generate_roots_from_simple(simple_roots, n_vars):
    """Generate full root system from simple roots via Weyl orbit BFS.

    simple_roots: list of tuples of Fractions, each of length n_vars.
    Returns list of all roots (as tuples), each unique.
    """
    def inner(a, b):
        return sum(x * y for x, y in zip(a, b))

    def reflect(beta, alpha):
        """Reflection of beta by alpha: beta - 2 (beta.alpha)/(alpha.alpha) alpha."""
        num = inner(beta, alpha)
        denom = inner(alpha, alpha)
        coef = 2 * num / denom
        return tuple(b - coef * a for b, a in zip(beta, alpha))

    # Seed with simple roots
    roots = set()
    frontier = list(simple_roots)
    for r in simple_roots:
        roots.add(r)
    while frontier:
        next_frontier = []
        for beta in frontier:
            for alpha in simple_roots:
                new_root = reflect(beta, alpha)
                if new_root not in roots:
                    roots.add(new_root)
                    next_frontier.append(new_root)
        frontier = next_frontier
    return list(roots)

def positive_roots(all_roots):
    """Extract positive roots: first non-zero coordinate is positive."""
    pos = []
    for r in all_roots:
        for x in r:
            if x != 0:
                if x > 0:
                    pos.append(r)
                break
    return pos

def e6_simple_roots():
    """E_6 simple roots in 6-dimensional basis.

    Using Bourbaki Planche V convention adapted to 6D:
      alpha_1 = (1, -1, 0, 0, 0, 0)
      alpha_2 = (0, 1, -1, 0, 0, 0)
      alpha_3 = (0, 0, 1, -1, 0, 0)
      alpha_4 = (0, 0, 0, 1, -1, 0)
      alpha_5 = (0, 0, 0, 0, 1, -1)
      alpha_6 = (1/2, 1/2, 1/2, 1/2, 1/2, 1/2) * sqrt(?) — but all roots have length sqrt(2),
               so we need alpha_6 chosen so that the Cartan relations hold.

    Simpler standard: use an 8D ambient and restrict.
    Let's use a 6D realization where E_6 is the root system in a subspace of
    R^8 with e_6 + e_7 + e_8 = 0, then project onto 6D.

    For computational simplicity: use A_5 + extra root in an auxiliary 6th dim.
      alpha_1 = (1, -1, 0, 0, 0, 0)
      alpha_2 = (0, 1, -1, 0, 0, 0)
      alpha_3 = (0, 0, 1, -1, 0, 0)
      alpha_4 = (0, 0, 0, 1, -1, 0)
      alpha_5 = (1/2, 1/2, 1/2, -1/2, -1/2, 1)       [short-ish root with 6th coord]
      alpha_6 = (1/2, 1/2, 1/2, 1/2, 1/2, 1)

    Check: all simple roots should have length sqrt(2). Length^2 of each:
      (1,-1,0,...): 2 ✓
      (1/2,...,1): 1/4+1/4+1/4+1/4+1/4+1 = 5/4+1 = 9/4 — NOT 2. Wrong.

    Let me use the correct 8D Bourbaki description and just keep all 8 coords.
    """
    # Bourbaki Planche V simple roots of E_6 in 8D:
    # alpha_1 = (1/2)(e_1 + e_8) - (1/2)(e_2 + e_3 + e_4 + e_5 + e_6 + e_7)
    # alpha_2 = e_1 + e_2
    # alpha_3 = e_2 - e_1
    # alpha_4 = e_3 - e_2
    # alpha_5 = e_4 - e_3
    # alpha_6 = e_5 - e_4
    zero = Fraction(0)
    half = Fraction(1, 2)
    mhalf = Fraction(-1, 2)
    alpha_1 = (half, mhalf, mhalf, mhalf, mhalf, mhalf, mhalf, half)
    alpha_2 = (Fraction(1), Fraction(1), zero, zero, zero, zero, zero, zero)
    alpha_3 = (Fraction(-1), Fraction(1), zero, zero, zero, zero, zero, zero)
    alpha_4 = (zero, Fraction(-1), Fraction(1), zero, zero, zero, zero, zero)
    alpha_5 = (zero, zero, Fraction(-1), Fraction(1), zero, zero, zero, zero)
    alpha_6 = (zero, zero, zero, Fraction(-1), Fraction(1), zero, zero, zero)
    return [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6]

def e7_simple_roots():
    """E_7 simple roots in 8D (Bourbaki Planche VI)."""
    zero = Fraction(0)
    half = Fraction(1, 2)
    mhalf = Fraction(-1, 2)
    alpha_1 = (half, mhalf, mhalf, mhalf, mhalf, mhalf, mhalf, half)
    alpha_2 = (Fraction(1), Fraction(1), zero, zero, zero, zero, zero, zero)
    alpha_3 = (Fraction(-1), Fraction(1), zero, zero, zero, zero, zero, zero)
    alpha_4 = (zero, Fraction(-1), Fraction(1), zero, zero, zero, zero, zero)
    alpha_5 = (zero, zero, Fraction(-1), Fraction(1), zero, zero, zero, zero)
    alpha_6 = (zero, zero, zero, Fraction(-1), Fraction(1), zero, zero, zero)
    alpha_7 = (zero, zero, zero, zero, Fraction(-1), Fraction(1), zero, zero)
    return [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7]

def compute_exceptional_c2_via_wick(G):
    """Uniform exact Wick-contraction c_2 computation for G in {E6, E7}.

    Simple roots in 8D; positive roots generated via Weyl orbit BFS.
    Polynomials in the full 8 variables; the <sum(alpha.w)^6>_mu term via Wick.
    """
    if G == 'E6':
        simples = e6_simple_roots()
        expected_count = 36
    elif G == 'E7':
        simples = e7_simple_roots()
        expected_count = 63
    else:
        raise ValueError(G)
    all_roots = generate_roots_from_simple(simples, 8)
    pos = positive_roots(all_roots)
    assert len(pos) == expected_count, f"Expected {expected_count} positive roots for {G}, got {len(pos)}"

    r, Rp, A_val, B_val, d, a, cp, lam = derived(G)
    assert Rp == expected_count

    # Consistency: A = sum |alpha|^2 over positive roots
    def norm_sq(alpha):
        return sum(x*x for x in alpha)
    A_check = sum(norm_sq(alpha) for alpha in pos)
    assert A_check == A_val, f"A mismatch for {G}: expected {A_val}, got {A_check}"

    # Build prod(alpha.w)^2 (polynomial in 8 variables)
    prod_poly = {tuple(0 for _ in range(8)): Fraction(1)}
    for alpha in pos:
        sq = linear_form_squared_poly(alpha)
        prod_poly = poly_multiply_nd(prod_poly, sq)

    # Build sum(alpha.w)^6
    barR_poly = {}
    for alpha in pos:
        p6 = linear_form_power_poly(alpha, 6)
        for k, c in p6.items():
            barR_poly[k] = barR_poly.get(k, Fraction(0)) + c
    barR_poly = {k: v for k, v in barR_poly.items() if v != 0}

    sigma_sq = Fraction(1, 2)
    Z = wick_int_nd(prod_poly, sigma_sq)
    num = wick_int_nd(poly_multiply_nd(prod_poly, barR_poly), sigma_sq)
    E_barR = num / Z

    # Gamma cross-checks
    def w_sq_poly_8():
        p = {}
        for i in range(8):
            key = tuple(2 if k == i else 0 for k in range(8))
            p[key] = Fraction(1)
        return p
    w_sq = w_sq_poly_8()
    # But caveat: in 8D, |w|^2 over the FULL 8 space is NOT the |w|^2 in the
    # r-dim root subspace. The Gaussian <|w|^4>_mu under
    # e^{-|w|^2 over 8D} weighted by prod(alpha.w)^2 is NOT the Gamma formula
    # -- because the Gaussian integrates over 8D, but the prod(alpha.w)^2
    # vanishes in directions orthogonal to the root span, so the integral
    # FACTORIZES into the root-span part (with rank r and measure mu) times
    # the complement part (Gaussian). For our observable <I_6> (which lives
    # in the root span), only the root-span Gaussian measure matters.
    #
    # In particular, <|w|^2>_mu (over the root span, unit-cov Gaussian) is a,
    # but <|w|^2>_{8D weighted} equals <w.w>_mu over root span + <complement>,
    # where <complement> = (8-r)*sigma_sq = (8-r)/2 (Gaussian in (8-r) dims).
    #
    # For simple-laced E_6/E_7/E_8 in 8D embedding, root span has dim r = 6,7,8.
    # For E_6: 2 extra dims. For E_7: 1 extra dim. For E_8: 0 extra dims.
    #
    # The extra Gaussian dims contribute a CONSTANT to <|w|^{2k}>_{8D weighted}.
    # The Wick integration automatically accounts for this (since we integrate
    # over 8D with the given weight).
    #
    # Our target quantities (<|w|^6>_{mu-over-root-span} etc.) SHOULD use the
    # SPAN-restricted measure, not the 8D weighted one. To get this, we project
    # the integrand onto the root span before integrating.
    #
    # For simplicity, we restrict simples to the root span (projecting to 6D/7D)
    # and compute in the intrinsic dimension r.
    # This requires computing the projection matrix, which is a non-trivial step.
    #
    # Alternative easier approach: directly compute <|w|^{2k}>_{mu} = a(a+1)...(a+k-1)
    # from the Gamma identity (which IS intrinsic to the root span).
    # We only need to exactly compute the <barR>_{mu} in the intrinsic measure.
    #
    # Claim: <sum(alpha.w)^6>_{mu-intrinsic} is still computable via 8D Wick,
    # because sum(alpha.w)^6 as a polynomial in the 8D w-variables is identical
    # whether we restrict to the root span or keep ambient 8D, AS LONG AS we
    # project |w|^{2k} terms appropriately.
    #
    # More carefully: the measure mu on the root span (dim r) has Gaussian
    # e^{-|w_r|^2} prod(alpha.w_r)^2 dw_r. Each alpha lies in the root span
    # (8D vector with some coords). For w_r in the root span embedded in 8D,
    # (alpha.w_r) agrees with (alpha.w_ambient) when w_ambient = w_r + 0 (in
    # orthogonal direction), which corresponds to our 8D Wick with sigma^2 = 1/2.
    # BUT the ambient 8D Gaussian also integrates over the 8-r orthogonal directions,
    # giving EXTRA factors (constants). These factor out as a CONSTANT.
    #
    # Specifically: Z_{8D} = Z_{root-span} * (pi)^{(8-r)/2} * (norm from sigma=1/2
    # in orthogonal dirs).
    # And <barR>_{8D} = <barR>_{root-span}, because barR is a polynomial purely
    # in the root span coords (alpha has 0 components orthogonal to root span,
    # so (alpha.w_ambient) = (alpha.w_span)). But Wick in 8D integrates over ALL
    # 8 coords, including orthogonal — those contribute only to Z and to
    # pure-coefficient factors.
    #
    # So: <barR>_{8D} / Z_{8D} = [<barR>_{root-span} * Z_{orth}] / [Z_{root-span} * Z_{orth}]
    #                          = <barR>_{root-span} / Z_{root-span} = <barR>_mu ✓
    #
    # CONCLUSION: our 8D Wick gives the correct <barR>_mu directly.
    # BUT cross-checks against Gamma moments require restricting to the root span,
    # which is more involved.  We skip those cross-checks for E_6/E_7.

    # Assemble c_2
    coef_barR = -Fraction(d*d, 360 * cp**3)
    coef_w6 = Fraction(5*d, 12 * (d+2)) - Fraction(1, 3) + Fraction((d+12)*d, 144*(d+2))
    coef_w4 = Fraction(d*d * (d+1), 288 * (d+2))
    coef_w8 = Fraction((d+12)**2, 288 * (d+2)**2)

    c2 = coef_barR * E_barR + coef_w6 * gamma_moment(a, 3) + coef_w4 * gamma_moment(a, 2) + coef_w8 * gamma_moment(a, 4)
    return c2, {'E_barR': E_barR, 'n_pos_roots': len(pos)}

def poly_multiply_nd(p, q):
    """Multiply two polynomials (dicts over tuple keys) in n variables."""
    result = {}
    for ka, ca in p.items():
        for kb, cb in q.items():
            key = tuple(a + b for a, b in zip(ka, kb))
            result[key] = result.get(key, Fraction(0)) + ca * cb
    return {k: v for k, v in result.items() if v != 0}

def wick_int_nd(poly, sigma_sq=Fraction(1, 2)):
    """Gaussian expectation in n dims, Cov(v_i, v_j) = sigma_sq * delta_ij.
    poly: dict {(k_1,...,k_n): coef}.
    """
    total = Fraction(0)
    for key, c in poly.items():
        if any(k % 2 for k in key):
            continue
        contribution = c
        for k in key:
            m = k // 2
            dd = Fraction(1)
            for i in range(1, m+1):
                dd *= (2*i - 1)
            contribution *= dd * sigma_sq**m
        total += contribution
    return total

def f4_positive_roots():
    """Return 24 positive roots of F_4 as tuples of Fractions in orthonormal basis."""
    roots = []
    # Long: +e_i + e_j and +e_i - e_j for 0 <= i < j <= 3, giving 12 roots
    for i in range(4):
        for j in range(i+1, 4):
            for s in (1, -1):
                v = [Fraction(0)] * 4
                v[i] = Fraction(1)
                v[j] = Fraction(s)
                roots.append(tuple(v))
    # Short: +e_i for i = 0..3, giving 4 roots
    for i in range(4):
        v = [Fraction(0)] * 4
        v[i] = Fraction(1)
        roots.append(tuple(v))
    # Short: (1/2)(e_0 + s_1 e_1 + s_2 e_2 + s_3 e_3) for s_i in ±1, giving 8 roots
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                roots.append((Fraction(1, 2), Fraction(s1, 2), Fraction(s2, 2), Fraction(s3, 2)))
    return roots

def linear_form_squared_poly(alpha):
    """(alpha . w)^2 as polynomial in (w_1, ..., w_n)."""
    n = len(alpha)
    poly = {}
    for i in range(n):
        for j in range(n):
            key = tuple((1 if k == i else 0) + (1 if k == j else 0) for k in range(n))
            c = alpha[i] * alpha[j]
            poly[key] = poly.get(key, Fraction(0)) + c
    return {k: v for k, v in poly.items() if v != 0}

def linear_form_power_poly(alpha, power):
    """(alpha . w)^power as polynomial in (w_1, ..., w_n)."""
    n = len(alpha)
    # Build linear form as poly with single-variable keys
    lin = {}
    for i in range(n):
        if alpha[i] != 0:
            key = tuple(1 if k == i else 0 for k in range(n))
            lin[key] = alpha[i]
    # Power via repeated multiplication
    result = {tuple(0 for _ in range(n)): Fraction(1)}
    for _ in range(power):
        result = poly_multiply_nd(result, lin)
    return result

def compute_F4_c2_via_wick():
    """Exact c_2(F_4) via symbolic Wick on full rank-6 invariant sum(alpha.w)^6.

    Formula:
      c_2 = -d^2/(360 c'^3) * <sum(alpha.w)^6>_mu
          + (5d/(12(d+2)) - 1/3) * <|w|^6>_mu
          + (d^2(d+1)/(288(d+2))) * <|w|^4>_mu
          + ((d+12)d/(144(d+2))) * <|w|^6>_mu
          + ((d+12)^2/(288(d+2)^2)) * <|w|^8>_mu

    where <...>_mu is under the weighted measure e^{-|w|^2} prod(alpha.w)^2 dw,
    which we compute as wick_int_nd(prod * target, 1/2) / wick_int_nd(prod, 1/2).

    For the Gamma-identity terms |w|^{2k}, we use the closed form a(a+1)...(a+k-1).
    For the <sum(alpha.w)^6>_mu term, we compute via Wick directly.
    """
    G = 'F4'
    r, Rp, A, B, d, a, cp, lam = derived(G)
    n_vars = r  # 4

    # Build prod(alpha.w)^2 as polynomial in (w_1, ..., w_4)
    roots = f4_positive_roots()
    assert len(roots) == Rp, f"Expected {Rp} roots, got {len(roots)}"

    prod_poly = {tuple(0 for _ in range(n_vars)): Fraction(1)}
    for alpha in roots:
        sq = linear_form_squared_poly(alpha)
        prod_poly = poly_multiply_nd(prod_poly, sq)

    # Build sum(alpha.w)^6 as polynomial
    barR_poly = {}
    for alpha in roots:
        p6 = linear_form_power_poly(alpha, 6)
        for k, c in p6.items():
            barR_poly[k] = barR_poly.get(k, Fraction(0)) + c
    barR_poly = {k: v for k, v in barR_poly.items() if v != 0}

    # Compute <prod> (normalizer) and <prod * barR>
    sigma_sq = Fraction(1, 2)
    Z = wick_int_nd(prod_poly, sigma_sq)
    num = wick_int_nd(poly_multiply_nd(prod_poly, barR_poly), sigma_sq)
    E_barR = num / Z

    # Cross-check: <prod * |w|^6> should match Gamma closed form
    # |w|^6 = (w_1^2 + w_2^2 + w_3^2 + w_4^2)^3
    def w_sq_poly(n):
        poly = {}
        for i in range(n):
            key = tuple(2 if k == i else 0 for k in range(n))
            poly[key] = Fraction(1)
        return poly
    w_sq = w_sq_poly(n_vars)
    w6 = poly_multiply_nd(poly_multiply_nd(w_sq, w_sq), w_sq)
    num_w6 = wick_int_nd(poly_multiply_nd(prod_poly, w6), sigma_sq)
    E_w6_wick = num_w6 / Z
    E_w6_gamma = gamma_moment(a, 3)
    assert E_w6_wick == E_w6_gamma, f"Gamma cross-check failed: {E_w6_wick} != {E_w6_gamma}"

    # Similarly for |w|^4, |w|^8
    w4 = poly_multiply_nd(w_sq, w_sq)
    w8 = poly_multiply_nd(w4, w4)
    E_w4_wick = wick_int_nd(poly_multiply_nd(prod_poly, w4), sigma_sq) / Z
    E_w8_wick = wick_int_nd(poly_multiply_nd(prod_poly, w8), sigma_sq) / Z
    assert E_w4_wick == gamma_moment(a, 2), "w4 cross-check failed"
    assert E_w8_wick == gamma_moment(a, 4), "w8 cross-check failed"

    # Now assemble c_2 using the master formula
    # c_2 = coef_barR * E_barR + coef_w6 * E_w6 + coef_w4 * E_w4 + coef_w8 * E_w8
    coef_barR = -Fraction(d*d, 360 * cp**3)
    coef_w6 = Fraction(5*d, 12 * (d+2)) - Fraction(1, 3) + Fraction((d+12)*d, 144*(d+2))
    coef_w4 = Fraction(d*d * (d+1), 288 * (d+2))
    coef_w8 = Fraction((d+12)**2, 288 * (d+2)**2)

    c2 = coef_barR * E_barR + coef_w6 * gamma_moment(a, 3) + coef_w4 * gamma_moment(a, 2) + coef_w8 * gamma_moment(a, 4)

    return c2, {
        'E_barR': E_barR,
        'coef_barR': coef_barR, 'coef_w6': coef_w6, 'coef_w4': coef_w4, 'coef_w8': coef_w8,
        'Ew4': gamma_moment(a, 2), 'Ew6': gamma_moment(a, 3), 'Ew8': gamma_moment(a, 4),
        'Z': Z, 'num_barR': num,
    }

def c2_assuming_I6_zero(G):
    """Compute c_2(G) exactly as a rational, ASSUMING <I_6(w)>_mu = 0.

    This reduces bar_R := sum(alpha.w)^6 to mu_1 * |w|^6 where
      mu_1 = B_3(G) * 15 / (r(r+2)(r+4))
    (from angular averaging: mu_1 is the radial coefficient).

    Then all expectations under mu are Gamma-factor closed-form.
    """
    from fractions import Fraction
    r, Rp, A, B, d, a, cp, lam = derived(G)
    # B_3 = sum |alpha|^6
    if G == 'G2':
        B_3 = Fraction(3 * 1 + 3 * 27)  # 3 short |alpha|^2=1, 3 long |alpha|^2=3
    elif G == 'F4':
        B_3 = Fraction(12 * 1 + 12 * 8)  # 12 short |^2=1, 12 long |^2=2
    else:  # E6, E7, E8 simply laced, all |alpha|^2 = 2
        B_3 = Fraction(Rp * 8)
    mu_1 = B_3 * Fraction(15, r * (r+2) * (r+4))
    # Assuming <I_6> = 0:
    # <bar R>_mu = mu_1 * <|w|^6>_mu = mu_1 * a(a+1)(a+2)
    Ew4 = gamma_moment(a, 2)
    Ew6 = gamma_moment(a, 3)
    Ew8 = gamma_moment(a, 4)
    # Formula coefficients (same as in compute_F4_c2_via_wick)
    coef_barR = -Fraction(d*d, 360 * cp**3)
    coef_w6 = Fraction(5*d, 12 * (d+2)) - Fraction(1, 3) + Fraction((d+12)*d, 144*(d+2))
    coef_w4 = Fraction(d*d * (d+1), 288 * (d+2))
    coef_w8 = Fraction((d+12)**2, 288 * (d+2)**2)
    c2_assumed = coef_barR * mu_1 * Ew6 + coef_w6 * Ew6 + coef_w4 * Ew4 + coef_w8 * Ew8
    return c2_assumed, mu_1, B_3

if __name__ == '__main__':
    # Verify c_1 closed form for each family
    print("=== c_1(G) verification against closed form -d(d+4)/16 ===")
    for G in ['G2', 'F4', 'E6', 'E7', 'E8']:
        c1, EP, EQ2 = c1_closed(G)
        r, Rp, A, B, d, a, cp, lam = derived(G)
        expected = Fraction(-d*(d+4), 16)
        match = "✓" if c1 == expected else "✗"
        print(f"  {G}: c_1 = {c1} (expected {expected}) {match}")

    # E_8: should give 346377539/45
    print("\n=== c_2(E_8) from rank-6 isotropy (no <I_6> term) ===")
    c2_E8, A6_E8, B6_E8 = c2_for_E8()
    expected_E8 = Fraction(346_377_539, 45)
    print(f"  c_2(E_8) computed = {c2_E8}")
    print(f"  c_2(E_8) expected = {expected_E8}")
    print(f"  A_6(E_8) = {A6_E8}, B_6(E_8) = {B6_E8}")
    match = "✓" if c2_E8 == expected_E8 else f"✗ (diff = {c2_E8 - expected_E8})"
    print(f"  match: {match}")
    print(f"  c_2(E_8) approx = {float(c2_E8):.4f}")

    # G_2: rigorous rational c_2
    print("\n=== c_2(G_2) via exact Wick on primary invariant I_6^{(0)} ===")
    c2_G2, diag = compute_G2_c2_via_wick()
    # Compare to paper's numerical 143.97
    print(f"  c_2(G_2) (exact rational) = {c2_G2} = {float(c2_G2):.6f}")
    print(f"  Paper's numerical value = 143.97")
    print(f"  Match to 3 sig figs: {abs(float(c2_G2) - 143.97) < 1.0}")
    print(f"  Diagnostics:")
    for k, v in diag.items():
        print(f"    {k} = {v}")

    # Discriminant check c_1^2 - 2 c_2 for G_2
    r, Rp, A, B, d, a, cp, lam = derived('G2')
    c1_G2 = Fraction(-d*(d+4), 16)
    disc = c1_G2**2 - 2 * c2_G2
    print(f"\n  c_1(G_2)^2 - 2 c_2(G_2) = {disc} = {float(disc):.4f}")
    print(f"  disc < 0 (series >= 1/2 for all n > 0): {disc < 0}")

    # ==================================================================
    # F_4: rank-4, 24 positive roots, exact Wick in 4 variables
    # ==================================================================
    print("\n=== c_2(F_4) via exact Wick on sum(alpha.w)^6 ===")
    c2_F4, diag_F4 = compute_F4_c2_via_wick()
    r, Rp, A, B, d, a, cp, lam = derived('F4')
    c1_F4 = Fraction(-d*(d+4), 16)
    disc_F4 = c1_F4**2 - 2 * c2_F4
    print(f"  c_2(F_4) (exact rational) = {c2_F4} = {float(c2_F4):.4f}")
    print(f"  Paper's numerical value = 17287.5")
    print(f"  c_1(F_4)^2 - 2 c_2(F_4) = {disc_F4} = {float(disc_F4):.4f}")
    print(f"  disc < 0 (series >= 1/2 for all n > 0): {disc_F4 < 0}")
    print(f"  <sum(alpha.w)^6>_mu = {diag_F4['E_barR']}")
    print(f"  coef_barR = {diag_F4['coef_barR']}")

    # ==================================================================
    # Closed-form c_2 ASSUMING <I_6>_mu = 0 (check per family)
    # ==================================================================
    print("\n=== Closed-form c_2(G) ASSUMING <I_6(w)>_mu = 0 ===")
    print(f"{'G':<4}{'mu_1':>15}{'B_3':>12}{'c_2_assumed':>30}{'vs ref':>15}")
    refs = {'G2': Fraction(165851, 1152), 'F4': Fraction(8401705, 486),
            'E6': 82867.0, 'E7': 680351.0, 'E8': Fraction(346377539, 45)}
    for G in ['G2', 'F4', 'E6', 'E7', 'E8']:
        c2_assumed, mu_1, B_3 = c2_assuming_I6_zero(G)
        ref = refs[G]
        ref_f = float(ref)
        c2_f = float(c2_assumed)
        diff_pct = 100 * (c2_f - ref_f) / ref_f
        print(f"{G:<4}{str(mu_1):>15}{str(B_3):>12}{str(c2_assumed):>30}{diff_pct:>12.3f}%")

    # ==================================================================
    # E_6 and E_7: same Wick procedure in principle, but the degree-72
    # and degree-126 polynomial products overwhelm direct symbolic
    # arithmetic.  A specialized Wick-contraction-over-pairings approach
    # (using inner products of the roots as a sparse graph) is required.
    # Left as a mechanical follow-up task.  Pass "--e6e7" to attempt
    # (will likely exhaust memory for non-trivial rank).
    # ==================================================================
    import sys as _sys
    if '--e6e7' in _sys.argv:
        for G_name, paper_val in [('E6', 82150.0), ('E7', 563000.0)]:
            print(f"\n=== c_2({G_name}) via exact Wick (Weyl orbit from simple roots in 8D) ===")
            try:
                c2_G, diag_G = compute_exceptional_c2_via_wick(G_name)
                r, Rp, A, B, d, a, cp, lam = derived(G_name)
                c1_G = Fraction(-d*(d+4), 16)
                disc_G = c1_G**2 - 2 * c2_G
                print(f"  positive roots found: {diag_G['n_pos_roots']} (expected {Rp})")
                print(f"  c_2({G_name}) (exact rational) = {c2_G} = {float(c2_G):.4f}")
                print(f"  Paper's numerical value = {paper_val}")
                print(f"  c_1({G_name})^2 - 2 c_2({G_name}) = {disc_G} = {float(disc_G):.4f}")
                print(f"  disc < 0 (series >= 1/2 for all n > 0): {disc_G < 0}")
            except Exception as e:
                print(f"  ERROR: {e}")
                import traceback
                traceback.print_exc()
