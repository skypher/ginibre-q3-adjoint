# FM-SEC65 (luna_max_mars) printed code: e = 1, C = 3, q = 0 family phi_2(h_(p+2) h_2^2 h_1^(p+2s+2)), normalized closed form for s >= 5.
from math import comb
from functools import lru_cache
import sympy as sp

# Direct Catalan-moment evaluator for the e=1 definition.
@lru_cache(None)
def moment(m, n):
    if m < n or (m-n) % 2:
        return 0
    j = (m-n)//2
    return (n+1)*comb(m, j)//((m+n)//2+1)

def row(a):
    return tuple(
        (comb(a, k) if 0 <= k <= a else 0)
        - (comb(a, k-1) if 0 <= k-1 <= a else 0)
        for k in range(a+2)
    )

def W(a, m, n):
    c = row(a)
    N = a+1
    return sum(c[k]*moment(N-k, m)*moment(k, n)
               for k in range(N+1))

def cg(m, n):
    return range(abs(m-n), m+n+1, 2)

def phi_direct(a, u, v, w):
    A, B, C = u+1, v+1, w+1
    val = 0
    for d in cg(A, B):
        for k in cg(d, C):
            val += W(a, k, 0)
        val -= W(a, d, C)
    for d in cg(A, C):
        val -= W(a, d, B)
    for d in cg(B, C):
        val -= W(a, d, A)
    return val

# Fast exact mirror-form evaluator for q=0, C=3.
def phi_q0(p, s):
    a = p+2*s+2
    N = a+1

    def low(k):
        return ((comb(a, k) if 0 <= k <= a else 0)
                - (comb(a, k-1) if 0 <= k-1 <= a else 0))

    def c(k):
        if k < 0 or k > N:
            return 0
        return -low(N-k) if 2*k > N else low(k)

    def D(k):
        return c(k)**2-c(k-1)*c(k+1)

    def P(x):
        return (sum(D(k) for k in range(x, x+4))
                - c(x)*c(x+3) + c(x-1)*c(x+4))

    def g(k):
        return c(k), c(k-1)

    def add(v, w):
        return v[0]+w[0], v[1]+w[1]

    def wedge(v, w):
        return v[0]*w[1]-v[1]*w[0]

    U = add(g(p+s), g(s))
    V = add(g(p+s+4), g(s-4))
    return P(p+s)-P(s-4)-wedge(U, V)

# Boundary values at s=3,4,5; compare both exact evaluators.
for s in (3, 4, 5):
    vals = []
    for p in range(11):
        a = p+2*s+2
        direct = phi_direct(a, p+2, 2, 2)
        assert direct == phi_q0(p, s)
        vals.append(direct)
    print("s=", s, "p=0..10 values=", vals,
          "minimum=", min(vals))

# First small outside-LS point and the two HT3 misses.
print("p=0,s=3 LS comparison:", 10*3, 7*0**2+51*0+202)
for a, p, s in ((24, 4, 9), (26, 14, 5)):
    direct = phi_direct(a, p+2, 2, 2)
    assert direct == phi_q0(p, s)
    print("witness (a,p,s), phi, 10s, LS RHS =",
          (a, p, s), direct, 10*s, 7*p*p+51*p+202)

# Symbolic normalized closed form for s>=5.
p, s = sp.symbols("p s", integer=True, nonnegative=True)
N = p+2*s+3

def r(j):
    ratio = sp.prod((p+s+9-m)/(s-5+m)
                    for m in range(1, j+6))
    return sp.factor((p+3-2*j)*ratio/N)

# c_(s+j) divided by binom(N,s-5).
def cl(j):
    return r(j)

# c_(p+s+t) reflected using c_(N-k)=-c_k.
def ch(t):
    return -cl(3-t)

def Dloc(f, k):
    return f(k)**2-f(k-1)*f(k+1)

def P_high():
    return (sum(Dloc(ch, k) for k in range(4))
            - ch(0)*ch(3)+ch(-1)*ch(4))

def P_low():
    f = lambda k: cl(k-4)
    return (sum(Dloc(f, k) for k in range(4))
            - f(0)*f(3)+f(-1)*f(4))

U = (ch(0)+cl(0), ch(-1)+cl(-1))
V = (ch(4)+cl(-4), ch(3)+cl(-5))
F = sp.factor(P_high()-P_low()-(U[0]*V[1]-U[1]*V[0]))

den = (s**2*(s-4)**2*(s-3)**2*(s-2)**2*(s-1)**2
       *(s+1)**2*(s+2)**2*(s+3)**2*(s+4)*N)
outer = (p+3)*(p+4)*(p+5)*(p+s+8)*(p+2*s+4)
Q = sp.cancel(F*den/outer)
poly = sp.Poly(Q, p, s)

assert poly.degree(p) == 12 and poly.degree(s) == 10
assert poly.coeff_monomial(p**12) == 1
assert poly.coeff_monomial(s**10) == 6000
assert poly.coeff_monomial(p**4*s**7) == -2800
print("normalized factor: phi/B^2 = outer*Q/den")
print("degrees of Q:", (poly.degree(p), poly.degree(s)))
print("Q coefficients (p^12, s^10, p^4*s^7):",
      (poly.coeff_monomial(p**12),
       poly.coeff_monomial(s**10),
       poly.coeff_monomial(p**4*s**7)))
