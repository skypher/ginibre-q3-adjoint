
from math import comb
from fractions import Fraction

def row(a, e):
    N = a + e
    c = [
        sum((-1)**i * comb(e, i) * comb(a, m-i)
            for i in range(e+1) if 0 <= m-i <= a)
        for m in range(N+1)
    ]
    return c, N

def det(p, q):
    return p[0]*q[1] - p[1]*q[0]

def sub(q, p):
    return (q[0]-p[0], q[1]-p[1])

cases = [
    (4, 56, 32, 10, 0),
    (4, 58, 33, 10, 0),
    (4, 59, 34, 10, 0),
    (4, 59, 33, 11, 1),
]

for r, a, x, C, edge in cases:
    e = 2*r - 3
    c, N = row(a, e)
    cf = lambda k: c[k] if 0 <= k <= N else 0
    g = [(cf(k), cf(k-1)) for k in range(x, x+C+2)]
    z = g[-1]

    assert all(det(g[i], g[i+1]) > 0 for i in range(len(g)-1))
    witness = next(i for i in range(1, len(g)-1)
                   if det(g[0], g[i]) < 0)
    assert det(g[0], z) > 0

    p, q = g[edge], g[edge+1]
    assert det(p, z) > 0 and det(z, q) > 0
    v = sub(q, p)
    denominator = det(z, v)
    t = Fraction(det(p, z), denominator)
    lam = Fraction(det(p, q), denominator)
    assert 0 < t < 1 and lam > 1
    assert (p[0] + t*v[0], p[1] + t*v[1]) == (lam*z[0], lam*z[1])

    D = [cf(k)**2 - cf(k-1)*cf(k+1) for k in range(N+1)]
    S = sum(D[x:x+C+1])
    T = det(g[0], z)
    print((r, a, x, C), 'edge=', edge,
          'negative-direction witness=', witness,
          'det(start,end)=', T,
          'det(start,witness)=', det(g[0], g[witness]))
    print('  t=', t, 'lambda=', lam, 'S-T=', S-T)
