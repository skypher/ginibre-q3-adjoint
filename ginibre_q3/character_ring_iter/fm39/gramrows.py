
import argparse
import sympy as S

parser = argparse.ArgumentParser(
    description="Derive symbolic Krawtchouk Gram rows and check their signs.")
parser.add_argument("--rows", type=int, default=4)
J = parser.parse_args().rows

N, t, i, z, h, u, v = S.symbols("N t i z h u v")
PN = lambda e: S.Poly(e, N, t, domain=S.QQ)
L = 4*J + 6

# Exact polynomial coefficients of the monic Krawtchouk polynomials.
kap = [S.Integer(1)]
for j in range(1, L//2):
    f = S.Poly(
        S.expand(i*(N-i+1)*kap[j-1].subs(t, i-1)), i, N)
    out = 0
    for (a, b), co in f.terms():
        power = (S.bernoulli(a+1, t)
                 - S.bernoulli(a+1, 1))/(a+1)
        out += co*N**b*power
    kap.append(S.expand(out))
    print("kappa", j, flush=True)

def choosepoly(q, j):
    return S.prod(q-a for a in range(j))/S.factorial(j)

f = [PN(0) for _ in range(L)]
g = [PN(0) for _ in range(L)]
for ell in range(L):
    if ell % 2 == 0:
        f[ell] = PN((-1)**(ell//2)*kap[ell//2])
    g[ell] = PN(sum(
        (-1)**j*kap[j]*2**(ell-2*j)
        *choosepoly(t-2*j, ell-2*j)
        for j in range(ell//2+1)))

def mul(a, b):
    return [
        sum((a[j]*b[d-j] for j in range(d+1)), PN(0))
        for d in range(L)]

def add(*arrays):
    return [
        sum((a[j] for a in arrays), PN(0))
        for j in range(L)]

def arr(expr):
    p = S.Poly(expr, z)
    return [PN(p.nth(j)) for j in range(L)]

na, ka = (N*z-1)/2, (N*z+1)/2
abar = add(
    mul(arr((na+z)*(ka+z)**2*(1+2*z)), mul(f, f)),
    mul(arr(na**2*(ka+2*z)), mul(g, g)),
    mul(arr(-(N-2*t)*z*(1+z)*na*(ka+z)), mul(f, g)))
assert all(abar[j].is_zero for j in range(3))

# Closed coefficient formula (1).
div = PN((N-t+1)*(N-t+2))
r = []
for ell in range(2*J+2):
    degree = 2*ell+3
    numerator = sum((
        abar[a]*PN(2**(degree-a)
                   *choosepoly(ell-t-1, degree-a))
        for a in range(degree+1)), PN(0))
    r.append(numerator.exquo(div))

cache = {}
def prodlead(row, offset, cross):
    key = row, offset, cross
    if key not in cache:
        d = t-2*row
        cache[key] = PN((-1)**offset*sum(
            kap[a].subs({N: N+2, t: d}, simultaneous=True)
            *kap[offset-a].subs(
                {N: N+2, t: d-2*cross}, simultaneous=True)
            for a in range(offset+1)))
    return cache[key]

aa, bb, A, B = [], [], [], []
for j in range(J+1):
    a = r[2*j]
    for q in range(j):
        a -= aa[q]*prodlead(q, 2*j-2*q, False)
        a -= 2*bb[q]*prodlead(q, 2*j-2*q-1, True)
    aa.append(a)

    b = r[2*j+1]
    for q in range(j+1):
        b -= aa[q]*prodlead(q, 2*j+1-2*q, False)
    for q in range(j):
        b -= 2*bb[q]*prodlead(q, 2*j-2*q, True)
    b = PN(b.as_expr()/2)
    bb.append(b)

    an = a.exquo(PN(S.prod(t-q for q in range(2*j))))
    bn = b.exquo(PN(S.prod(t-q for q in range(2*j+2))))
    A.append(S.expand(an.as_expr().subs(N, h+t-2)))
    B.append(S.expand(bn.as_expr().subs(N, h+t-2)))

    print("ROW", j, "A =", S.collect(A[j], h), flush=True)
    print("ROW", j, "B =", S.collect(B[j], h), flush=True)

for j in range(1, J+1):
    F = (h*A[j] + (t-2*j+2)*h*h*B[j-1]
         + (t-2*j-1)*B[j])
    tests = [
        (-B[j]).subs(
            {t: 2*j+2+u, h: 2*j+4+u+v}, simultaneous=True),
        F.subs(
            {t: 2*j+1+u, h: 2*j+3+u+v}, simultaneous=True),
        (A[j]+2*h*B[j-1]).subs(
            {t: 2*j, h: 2*j+2+v}, simultaneous=True)]
    stats = []
    for f in tests:
        p = S.Poly(f.expand(), u, v)
        assert all(co > 0 for co in p.coeffs())
        stats.append((len(p.terms()), min(p.coeffs())))
    print("row", j, stats, flush=True)

# Obstruction to the specified two-step Gram transfer.
if J >= 2:
    remainder = S.expand(
        A[2] - h*(h-1)*A[1].subs(
            {t: t-2, h: h+2}, simultaneous=True))
    boundary = S.factor(remainder.subs(h, t+2))
    print("two-step remainder:", boundary)
    assert boundary.subs(t, 64) == -S.Rational(583692, 7)
