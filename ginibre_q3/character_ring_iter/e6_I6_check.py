#!/usr/bin/env python3
"""Numerical check of <I_6(E_6)>_mu = 0 via outer-automorphism argument.

The E_6 outer automorphism sigma is the nontrivial element of
Out(E_6) = Z/2.  Induced by the Dynkin-diagram involution swapping
nodes alpha_1 <-> alpha_6, alpha_3 <-> alpha_5, alpha_2, alpha_4 fixed.

As an orthogonal transformation on the Cartan (8D Bourbaki embedding of E_6):
sigma acts specifically.  We compute sigma as a matrix, verify it permutes
the root set, and check that the natural primary I_6 candidate is
sigma-anti-invariant.

A CANDIDATE primary I_6 for E_6:  from the invariant theory, I_6 is the
degree-6 polynomial that is NOT expressible via I_2^3 alone.  One
construction: compute the adjoint trace
  I_6^{can}(w) := tr((H(w))^6) - c * (tr(H(w)^2))^3
where H(w) = sum_i w_i h_i (the Cartan element), and c is chosen so
that the tr(H^2)^3 subtraction removes the I_2^3 component.  For E_6
in 8D embedding, an explicit construction uses the 27-dim fundamental
rep.
"""
import numpy as np

def e6_simples():
    zero = 0.0
    half = 0.5
    mhalf = -0.5
    return np.array([
        [half, mhalf, mhalf, mhalf, mhalf, mhalf, mhalf, half],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
    ], dtype=np.float64)

def construct_sigma():
    """E_6 outer automorphism sigma as an orthogonal 8x8 matrix.

    sigma acts on simple roots: alpha_1 <-> alpha_6, alpha_2 fixed,
    alpha_3 <-> alpha_5, alpha_4 fixed.
    """
    simples = e6_simples()
    # Target images under sigma: alpha_1 <-> alpha_6, alpha_3 <-> alpha_5, alpha_2, alpha_4 fixed
    # Assuming Bourbaki E_6 ordering.
    image_idx = [5, 1, 4, 3, 2, 0]  # alpha_1 -> alpha_6, alpha_2 -> alpha_2, etc.
    # sigma as linear map: sigma(simples[i]) = simples[image_idx[i]]
    # Solve for sigma: sigma * S^T = S_img^T (where S is matrix of simples by rows)
    S = simples  # (6, 8)
    S_img = simples[image_idx]  # (6, 8)
    # Want sigma such that sigma @ S.T = S_img.T, i.e., S @ sigma.T = S_img, i.e., sigma = S_img.T @ pinv(S.T)
    # Pseudoinverse since sigma is 8x8 but simples is 6x8 (6D in 8D)
    # sigma restricted to root span is determined.  On orthogonal complement (2D in 8D),
    # sigma can be ANY orthogonal transformation.  Pick identity on orthogonal complement.
    # Use least-squares: sigma = S_img.T @ pinv(S.T) is 8x8 but only constrains 6D subspace.
    # More careful: sigma * S_i = S_img_i for each i. So S @ sigma.T = S_img, giving
    # sigma.T = pinv(S) @ S_img (min-norm solution).  On orthogonal complement sigma.T projects to 0.
    # Then add identity on orthogonal complement.
    u, s, vt = np.linalg.svd(S, full_matrices=True)
    # Root span is spanned by first 6 rows of V (vt[:6].T basis).
    # Orthogonal: vt[6:].T
    V_span = vt[:6]  # (6, 8)
    V_perp = vt[6:]  # (2, 8)
    # Project simples to span basis: S_span_coords = S @ V_span.T  (6, 6)
    S_span = S @ V_span.T  # (6, 6)
    S_img_span = S_img @ V_span.T  # (6, 6)
    # sigma_span: S_span @ sigma_span.T = S_img_span, so sigma_span.T = inv(S_span) @ S_img_span
    sigma_span = np.linalg.solve(S_span, S_img_span).T  # Check: S_span @ sigma_span.T = S_img_span
    # Now build full sigma as acting sigma_span on V_span coords, identity on V_perp coords
    sigma = V_span.T @ sigma_span @ V_span + V_perp.T @ V_perp
    return sigma

def verify_sigma_orthogonal(sigma):
    """sigma should be an orthogonal matrix (preserves inner products)."""
    err = np.linalg.norm(sigma.T @ sigma - np.eye(8))
    return err

def generate_positive_roots():
    """Generate E_6 positive roots via Weyl orbit BFS."""
    simples = e6_simples()
    tol = 1e-9
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
                nv = b - num * alpha
                k = key(nv)
                if k not in roots:
                    roots[k] = nv.tolist()
                    new_frontier.append(nv.tolist())
        frontier = new_frontier
    all_r = np.array(list(roots.values()))
    # Positive: first nonzero coord positive
    out = []
    for r in all_r:
        for x in r:
            if abs(x) > tol:
                if x > 0:
                    out.append(r)
                break
    return np.array(out), all_r

def verify_sigma_preserves_roots(sigma, all_roots):
    """Check sigma(alpha) is a root for every root alpha."""
    tol = 1e-7
    ok = True
    root_set = {tuple(round(float(x), 6) for x in r) for r in all_roots}
    for r in all_roots:
        r_img = sigma @ r
        key = tuple(round(float(x), 6) for x in r_img)
        if key not in root_set:
            ok = False
            print(f"  Root {r} maps to {r_img}, not in root set!")
            break
    return ok

def construct_I6_candidate(pos_roots):
    """A candidate primary degree-6 E_6 invariant.

    I_6(w) := sum_{alpha in R+}(alpha.w)^6 - mu_1 * |w|^6
    where mu_1 = angular average = B_3 * 15 / (r(r+2)(r+4)).

    This is the `radial-subtracted` part of sum(alpha.w)^6.
    It's Weyl-invariant (as the difference of two W-invariants).
    """
    r, Rp = 6, 36  # rank, |R_+|
    # B_3 = sum |alpha|^6; for simply-laced |alpha|^2 = 2: B_3 = 36 * 8 = 288
    B_3 = 36 * 8
    mu_1 = B_3 * 15 / (r * (r+2) * (r+4))
    def I6(w):
        bar_R = np.sum(np.sum(pos_roots * w, axis=1)**6)
        w_sq_3 = np.sum(w * w)**3
        return bar_R - mu_1 * w_sq_3
    return I6, mu_1

def main():
    simples = e6_simples()
    pos_roots, all_roots = generate_positive_roots()
    print(f"E_6 positive roots: {len(pos_roots)}, total: {len(all_roots)}")
    assert len(pos_roots) == 36

    sigma = construct_sigma()
    ortho_err = verify_sigma_orthogonal(sigma)
    print(f"sigma orthogonality error: {ortho_err:.2e}")

    preserves = verify_sigma_preserves_roots(sigma, all_roots)
    print(f"sigma preserves root set: {preserves}")

    I6, mu_1 = construct_I6_candidate(pos_roots)
    print(f"\nmu_1 = {mu_1}")

    # Check: is I_6 sigma-invariant or anti-invariant?
    # Sample random points and compare I_6(w) with I_6(sigma @ w).
    print("\nChecking sigma-action on I_6 candidate bar_R - mu_1 |w|^6:")
    rng = np.random.default_rng(0)
    for trial in range(5):
        w = rng.standard_normal(8)
        w_sigma = sigma @ w
        I_val = I6(w)
        I_sigma = I6(w_sigma)
        ratio = I_sigma / I_val if abs(I_val) > 1e-12 else 0
        print(f"  w = {w[:3]}..., I_6(w) = {I_val:.4f}, I_6(sigma w) = {I_sigma:.4f}, ratio = {ratio:.4f}")

    # If I_6 is sigma-invariant (ratio = 1), outer aut doesn't force <I_6> = 0.
    # If anti-invariant (ratio = -1), <I_6>_mu = 0 by sigma-invariance of mu.
    # If some other relation, further analysis needed.

    # Alternative: check if the MC-estimated <I_6>_mu is consistent with 0.
    print("\nSampling <I_6(w)>_mu via importance sampling:")
    import math
    N = 2_000_000
    batch = 100_000
    total_num = 0.0
    total_denom = 0.0
    R = pos_roots
    for _ in range(N // batch):
        W = rng.standard_normal(size=(batch, 8))
        inner = W @ R.T
        weights = np.prod(inner * inner, axis=1) * np.exp(-(W * W).sum(axis=1) / 2)
        bar_R = np.sum(inner**6, axis=1)
        w_sq = (W * W).sum(axis=1)
        w_sq_3 = w_sq ** 3
        I6_vals = bar_R - mu_1 * w_sq_3
        total_num += (I6_vals * weights).sum()
        total_denom += weights.sum()
    mean_I6 = total_num / total_denom
    # Stderr rough
    print(f"  <I_6(E_6)>_mu (MC, N={N}) ≈ {mean_I6:.4f}")
    # Compare to <|w|^6>_mu = a(a+1)(a+2) = 39*40*41 = 63960, so mu_1*<|w|^6> = 9*63960 = 575640
    # Scale of bar_R is similar
    scale = mu_1 * 63960
    print(f"  Relative to mu_1 * <|w|^6> = {scale}: {mean_I6 / scale * 100:.3f}%")

if __name__ == '__main__':
    main()
