# FM-SEC45 (luna_max_uranus) gamma <= N certificate, verbatim helper functions (the check() asserts replaced in gammaN_cert_census.py)
from math import comb, isqrt
from fractions import Fraction

def c_row(a, e):
    return [sum((-1)**h * comb(e, h) * comb(a, k-h)
                for h in range(e+1) if 0 <= k-h <= a)
            for k in range(a+e+1)]

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def g(c, k):
    return at(c, k), at(c, k-1)

def add(u, v):
    return u[0] + v[0], u[1] + v[1]

def wedge(u, v):
    return u[0]*v[1] - u[1]*v[0]

def D(c, k):
    return at(c, k)**2 - at(c, k-1)*at(c, k+1)

def B(c, k):
    return at(c, k-1) + at(c, k+1)

def Ac(c, k):
    return at(c, k-1) - at(c, k+1)

def P(c, x, C):
    T = at(c, x)*at(c, x+C) - at(c, x-1)*at(c, x+C+1)
    return sum(D(c, k) for k in range(x, x+C+1)) - T

def upoly(n):
    if n == 0:
        return [1]
    if n == 1:
        return [0, 1]
    u0, u1 = [1], [0, 1]
    for _ in range(2, n+1):
        u2 = [0]*(len(u1)+1)
        for k, v in enumerate(u1):
            u2[k+1] += v
        for k, v in enumerate(u0):
            u2[k] -= v
        while u2 and u2[-1] == 0:
            u2.pop()
        u0, u1 = u1, u2
    return u1

def mul(f, h):
    out = {}
    for (i, j), x in f.items():
        for (k, l), y in h.items():
            out[i+k, j+l] = out.get((i+k, j+l), 0) + x*y
    return {ij: v for ij, v in out.items() if v}

def phi_direct(e, a, A, B, C):
    f = {(e-k, k): (-1)**k * comb(e, k) for k in range(e+1)}
    f = mul(f, {(a-k, k): comb(a, k) for k in range(a+1)})
    for n in (A, B, C):
        h = {}
        for k, v in enumerate(upoly(n)):
            h[k, 0] = h.get((k, 0), 0) + v
            h[0, k] = h.get((0, k), 0) - v
        f = mul(f, {ij: v for ij, v in h.items() if v})

    def moment(m):
        return 0 if m % 2 else comb(m, m//2) // (m//2 + 1)

    z = sum(v*moment(i)*moment(j) for (i, j), v in f.items())
    assert z % 2 == 0
    return z // 2

def ceil_sqrt(x):
    n = isqrt(x.numerator // x.denominator)
    return n if n*n*x.denominator == x.numerator else n+1

def check(e, a, C, p, q, s):
    N = a + e
    A, B0 = C+p, C+q
    V0 = (a+1)*(e+1)
    c = c_row(a, e)
    assert e % 2 == 1 and C >= 3 and p >= q >= 0 and s >= 0
    assert N == C+p+q+2*s

    alpha, tau = p+s, N-s+1
    U = add(g(c, alpha), g(c, q+s))
    V = add(g(c, tau), g(c, s-C-1))
    X = wedge(U, V)
    delta = P(c, alpha, C) - P(c, s-C-1, C)

    terms = []
    Lambda = 0
    for d in range(abs(A-B0), A+B0+1, 2):
        j0 = (N+d-C)//2
        i = j0+C+1
        if d < C:
            j, eps = N-j0, 1
        else:
            j, eps = j0, -1

        E = D(c, j)-D(c, i)
        W = at(c, j)*B(c, i)-B(c, j)*at(c, i)
        candidates = []
        metric_R, root = None, 0

        if eps*W >= 0:
            lam, source = E, "OL"
        else:
            x, y = 2*j-N, 2*i-N
            if 0 < y*y < 4*V0:
                rhs = D(c, i)*(4*y*y*D(c, j)
                               -(y*y-x*x)*Ac(c, j)**2)
                if rhs >= 0:
                    R = Fraction(rhs, 4*V0-y*y)
                    if E*E >= R:
                        candidates.append((R, "Mi"))
            if 0 < x*x < 4*V0:
                rhs = D(c, j)*(4*x*x*D(c, i)
                               +(y*y-x*x)*Ac(c, i)**2)
                if rhs >= 0:
                    R = Fraction(rhs, 4*V0-x*x)
                    if E*E >= R:
                        candidates.append((R, "Mj"))

            assert candidates
            metric_R, source = min(candidates)
            assert W*W <= metric_R
            root = ceil_sqrt(metric_R)
            lam = E-root
            assert lam >= 0

        step = E+eps*W
        assert step >= lam
        terms.append((d, j, i, E, W, eps, step, lam,
                      source, metric_R, root))
        Lambda += lam

    assert sum(t[6] for t in terms) == delta
    assert phi_direct(e, a, A, B0, C) == delta-X

    d0, sigma = a-e, N+2
    au = 4*(U[0]**2-U[1]**2)
    bu = 4*(d0*U[0]-sigma*U[1])**2
    av = 4*(V[0]**2-V[1]**2)
    bv = 4*(d0*V[0]-sigma*V[1])**2
    K = 16*Lambda**2 + au*av
    L = 64*V0*Lambda**2 - au*bv - bu*av
    Z = bu*bv

    assert K > 0 and 0 < L < 8*V0*K and L*L >= 4*K*Z
    t = Fraction(L, 2*K)
    margin = Fraction(L*L-4*K*Z, 4*K)
    assert 0 < t < 4*V0 and margin >= 0

    return ((e, a, C, p, q, s), (A, B0, C), delta, X, delta-X,
            Lambda, (t, 4*V0, margin), terms)

