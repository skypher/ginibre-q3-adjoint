
import argparse
import sympy as S

parser = argparse.ArgumentParser()
parser.add_argument("--rows", type=int, default=7)
args = parser.parse_args()

x, N, v = S.symbols("x N v")
domain = S.QQ.poly_ring(N)

def Q(e):
    return S.Poly(e, x, domain=domain)

def binomial_poly(a, b):
    if b < 0:
        return S.Integer(0)
    return S.prod(a-i for i in range(b))/S.factorial(b)

def kraw(M, degree):
    K = [Q(1), Q(x)]
    for d in range(1, degree):
        K.append(Q(x)*K[d]-d*(M-d+1)*K[d-1])
    return K

def direct(t):
    old = kraw(N, t)
    K = kraw(N+2, 2*t)
    f, g = old[t], old[t].shift(2)
    n, k = Q((N-x)/2), Q((N+x)/2)
    numerator = (
        (n+1)*(k+1)**2*Q(x+2)*f**2
        + Q(x)*n**2*(k+2)*g**2
        - (N-2*t)*Q(x+1)*n*(k+1)*f*g
    )
    rem = numerator.exquo(Q((x+1)*(N-t+1)*(N-t+2)))
    rr = [S.Integer(0)]*(t+1)
    for d in range(t, -1, -1):
        rr[d] = rem.nth(2*d)
        rem -= Q(rr[d])*Q(x*(x+2))**d
    assert rem.is_zero
    R = Q(sum(rr[d]*x**(2*d) for d in range(t+1)))
    rem, aa, bb = R, [], []
    for j in range(t//2+1):
        aa.append(S.expand(rem.nth(2*t-4*j)))
        rem -= Q(aa[-1])*K[t-2*j]**2
        if j < t//2:
            bb.append(S.expand(rem.nth(2*t-4*j-2)/2))
            rem -= 2*Q(bb[-1])*K[t-2*j]*K[t-2*j-2]
    assert rem.is_zero
    return R, K, aa, bb

def source_coefficients(R, K, t, M):
    rem = R
    g = [S.Integer(0)]*(t+1)
    for r in range(t, -1, -1):
        g[r] = rem.nth(2*r)
        rem -= Q(g[r])*K[2*r]
    assert rem.is_zero

    c = []
    for d in range(t+1):
        value = S.Integer(0)
        for r in range(d, t+1):
            exponent = 2*r-M-1
            kernel = (
                binomial_poly(exponent, r-d)
                - binomial_poly(exponent, r-d-1)
            )
            value += S.factorial(r)**2*g[r]*kernel
        c.append(S.expand(value/S.factorial(d)**2))

    rem = R
    for d in range(t, -1, -1):
        assert S.expand(rem.nth(2*d)-c[d]) == 0
        rem -= Q(c[d])*K[d]**2
    assert rem.is_zero
    return c

def gram_from_source(c, t, M):
    H = S.Integer(0)
    aa, bb = [], []
    for j in range(t//2+1):
        d = t-2*j
        aa.append(S.expand(c[d]+(M-2*d)*H))
        if j < t//2:
            a = d*(M-d+1)
            bb.append(S.expand((c[d-1]+a*(M-2*d+2)*H)/2))
            H = S.expand(c[d-1]+a*a*H)
    return aa, bb

def terminal_from_source(c, t, M, h):
    total = c[0]
    for ell in range(t//2):
        if ell == 0:
            weight = M+h
        else:
            weight = (
                2*(M-1)*(2*M*(M-1)+h*(M-2))
                * S.prod(
                    (2*q*(M-2*q+1))**2
                    for q in range(2, ell+1)
                )
            )
        total += c[2*ell+1]*weight
    return S.expand(total/S.factorial(t))

source_checks = gram_checks = 0
levels = list(range(2, 2*args.rows+1, 2))
levels += [2*args.rows+2, 2*args.rows+3]

for t in levels:
    M = N+2
    R, K, aa, bb = direct(t)
    c = source_coefficients(R, K, t, M)
    anew, bnew = gram_from_source(c, t, M)
    assert all(S.expand(a-b) == 0 for a,b in zip(aa,anew))
    assert all(S.expand(a-b) == 0 for a,b in zip(bb,bnew))
    source_checks += len(c)
    gram_checks += len(aa)+len(bb)

    if t <= 2*args.rows:
        j, h = t//2, N-t+2
        z = terminal_from_source(c, t, M, h)
        expected = (aa[j]+2*h*bb[j-1])/S.factorial(t)
        assert S.expand(z-expected) == 0
        assert S.Poly(z,N).rem(S.Poly(N+2,N)).is_zero
        zv = S.Poly(z.subs(N,2*t+v),v)
        print("Z",j,"terms",len(zv.terms()),
              "min",min(zv.coeffs()),flush=True)
    else:
        print("both source and Gram formulas agree at t =",t,
              flush=True)

    if t == 8:
        print("source c_1 at t=8:",
              S.expand(c[1].subs(N,16+v)),flush=True)

print("source coefficient checks:",source_checks)
print("Gram coefficient checks:",gram_checks)
