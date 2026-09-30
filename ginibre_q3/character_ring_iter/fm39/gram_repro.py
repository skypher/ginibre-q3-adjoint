
import sympy as S
from math import comb, factorial
from fractions import Fraction as F

p, N, T, X, v = S.symbols("p N T X v")
P = lambda e: S.Poly(e, p, N, domain=S.QQ)
Q = lambda e: S.Poly(e, X, N, domain=S.QQ)
PN = lambda e: S.Poly(e, N, domain=S.QQ)

n, k = P((N-p)/2), P((N+p)/2)
K, B = [P(1), P(p)], [Q(1), Q(X)]
for j in range(1, 16):
    K.append(P(p)*K[-1] - P(j*(N-j+1))*K[-2])
    B.append(Q(X)*B[-1] - Q(j*(N+3-j))*B[-2])

def Rpoly(t):
    f = K[t]
    g = P(f.as_expr().subs(p, p+2).expand())
    A = ((n+1)*(k+1)**2*(p+2)*f*f
         + p*n*n*(k+2)*g*g
         - (N-2*t)*(p+1)*n*(k+1)*f*g)
    A = A.exquo(P((p+1)*(N-t+1)*(N-t+2)))
    A = S.Poly(A.as_expr(), p, N, T).rem(
        S.Poly(p*p+2*p-T, p, N, T))
    assert A.degree(p) == 0
    return S.Poly(A.as_expr(), T, N)

def cx(f, d):
    return S.Poly.from_dict(
        {(j,): c for (i, j), c in f.terms() if i == d},
        N, domain=S.QQ)

def gram(t, R):
    rem = Q(R.as_expr().subs(T, X*X))
    aa, bb = [], []
    for d in range(t, t % 2 - 1, -2):
        a = cx(rem, 2*d)
        aa.append(a)
        rem -= Q(a.as_expr())*B[d]**2
        if d >= 2:
            b = PN(cx(rem, 2*d-2).as_expr()/2)
            bb.append(b)
            rem -= Q(2*b.as_expr())*B[d]*B[d-2]
    assert rem.is_zero
    return aa, bb

C = lambda n, k: comb(n, k) if 0 <= k <= n else 0

def c(N, t, k):
    return sum((-1)**j*C(t, j)*C(N-t, k-j)
               for j in range(t+1))

direct = signpolys = 0
for t in range(2, 17):
    R = Rpoly(t)
    aa, bb = gram(t, R)
    h = PN(N-t+2)

    assert aa[0] == PN(1) and bb[0].is_zero
    a1 = t*(t-1)*(
        (N-t+2)*(N-t+3) - S.Rational(2, 3)*(t+1))
    assert aa[1] == PN(a1)
    if len(bb) > 1:
        b1 = (-S.Rational(2, 15)*t*(t-1)*(t-2)*(t-3)
              *(5*N-2*t+18))
        assert bb[1] == PN(b1)

    targets = [-b for b in bb]
    for j, a in enumerate(aa):
        row = a
        if j:
            row += bb[j-1]*(t-2*j+2)*h
        if j < len(bb):
            row = row*(t-2*j)*h + bb[j]
        targets.append(row)

    for f in targets:
        shifted = S.Poly(
            f.as_expr().subs(N, 2*t+v).expand(), v)
        assert all(co >= 0 for co in shifted.coeffs())
        signpolys += 1

    if t in (4, 5, 6, 7, 8, 10):
        for NN in range(t, t+13):
            fall = factorial(NN)//factorial(NN-t)
            for kk in range(NN//2+1, NN+1):
                pp, nn = 2*kk-NN, NN-kk
                z = [c(NN, t, kk+j) for j in (-1, 0, 1, 2)]
                actual = (z[1]**2-z[0]*z[2]
                          -z[2]**2+z[1]*z[3])
                got = F(
                    C(NN, kk)**2*(pp+1)*(NN-t+1)*(NN-t+2)
                    *int(R.eval({N: NN, T: pp*(pp+2)})),
                    (nn+1)*(kk+1)**2*(kk+2)*fall**2)
                assert actual == got
                direct += 1

t = S.symbols("t", integer=True)
h = N-t+2
a1 = t*(t-1)*(h*(h+1)-S.Rational(2, 3)*(t+1))
b1 = (-S.Rational(2, 15)*t*(t-1)*(t-2)*(t-3)
      *(5*N-2*t+18))
E = S.cancel(((t-2)*h*a1+b1)/(t*(t-1)*(t-2)))
E = S.Poly(E.subs(N, 2*t+v).expand(), t, v)
assert all(co > 0 for co in E.coeffs())

print("Direct checks:", direct)
print("Coefficientwise sign-polynomial checks:", signpolys)
print("Symbolic row-one polynomial:", E.as_expr())
