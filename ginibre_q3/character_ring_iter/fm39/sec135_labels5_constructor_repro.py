from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from datetime import datetime, timezone
import sympy as sp

q, z = sp.symbols("q z")
rho, delta = Q(6, 7), Q(1, 7)
tau, J = Q(17, 21), Q(352, 1225)

def bernstein_matrix(deg, lo, hi):
    return [[
        sum(Q(comb(i,t))*lo**(i-t)*(hi-lo)**t
            * Q(comb(k,t),comb(deg,t))
            for t in range(min(i,k)+1)
        )
        for i in range(deg+1)]
        for k in range(deg+1)]

def cert(expr, zlo=Q(0), zhi=rho, nq=16, nz=16):
    poly = sp.Poly(sp.expand(expr), q, z)
    dq, dz = poly.degree(q), poly.degree(z)
    a = [[Q(poly.coeff_monomial(q**i*z**j))
          for j in range(dz+1)] for i in range(dq+1)]
    tz = [bernstein_matrix(
              dz, zlo+(zhi-zlo)*j/nz, zlo+(zhi-zlo)*(j+1)/nz)
          for j in range(nz)]
    minimum, negative_boxes = None, 0
    for iq in range(nq):
        tq = bernstein_matrix(
            dq, Q(-1)+Q(2*iq,nq), Q(-1)+Q(2*(iq+1),nq))
        part = [[sum(tq[k][i]*a[i][j] for i in range(dq+1))
                 for j in range(dz+1)] for k in range(dq+1)]
        for T in tz:
            vals = [sum(part[k][j]*T[l][j] for j in range(dz+1))
                    for k in range(dq+1) for l in range(dz+1)]
            bmin = min(vals)
            minimum = bmin if minimum is None else min(minimum,bmin)
            negative_boxes += (bmin < 0)
    return minimum, negative_boxes

u = 7*(z-Q(6,7))
F3p = 4*z*(1+q+q**2)-2
F3m = 4*z*(1-q+q**2)-2
for name,F in (("Z+P",F3p),("Z-P",F3m)):
    D = sp.expand(F.subs(z,1))
    checks = {
        "weighted": (Q(6,7)**3)*D**2-z*F**2,
        "raw lower": 2*D+F,
        "raw upper": 2*D-F,
        "positive u": F-D*u,
    }
    for label,expr in checks.items():
        interval = (rho, Q(1)) if label == "positive u" else (Q(0),rho)
        assert cert(expr,*interval)[1] == 0

F4m = 4*(1+q**2)*z-3
D4m = sp.expand(F4m.subs(z,1))
for label,expr,lo,hi in (
    ("weighted",Q(6,7)**4*D4m**2-z**2*F4m**2,Q(0),rho),
    ("raw lower",3*D4m+F4m,Q(0),rho),
    ("raw upper",3*D4m-F4m,Q(0),rho),
    ("positive u",F4m-D4m*u,rho,Q(1)),
):
    assert cert(expr,lo,hi)[1] == 0

F4p = 16*(1+q**4)*z**2-12*(1+q**2)*z+2
D4p = sp.expand(F4p.subs(z,1))
theta4 = Q(85,147)
checks4 = (
    (theta4*D4p+F4p,Q(0),Q(7,9)),
    (theta4*D4p-F4p,Q(0),Q(7,9)),
    (F4p,Q(7,9),rho),
    (theta4*D4p-F4p,Q(7,9),rho),
    (F4p-D4p*u**2,rho,Q(1)),
)
for expr,lo,hi in checks4:
    assert cert(expr,lo,hi)[1] == 0

A = 1+q+q**2+q**3+q**4
B = 1+q+q**2
F5 = 16*A*z**2-16*B*z+3
D5 = sp.expand(F5.subs(z,1))
assert sp.expand(
    D5-Q(21,16)-(4*q+3)**2*(16*q**2-8*q+3)/16
) == 0
checks5 = (
    (F5.subs(z,Q(6,7)),Q(0),rho),
    (2*Q(6,7)*16*A-16*B,Q(0),Q(1)),
    (Q(9,16)*D5**2-z*F5**2,Q(0),rho),
    (3*D5+F5,Q(0),rho),
    (3*D5-F5,Q(0),rho),
    (F5-D5*u**2,rho,Q(1)),
)
for expr,lo,hi in checks5:
    assert cert(expr,lo,hi)[1] == 0
# q -> -q gives the same bounds for W5(-q).

C = Q(36369,200)
def tail(h):
    return C*(Q(15*h,7)+1)*(Q(15*h,7)+2)

b68 = tail(68)*rho**99
b67_squared = tail(67)**2*rho**195
ratio67_squared = (
    rho**3*Q(15*67+22,15*67+7)**2
    * Q(15*67+29,15*67+14)**2
)
assert b68 < 1 and b67_squared > 1 and ratio67_squared < 1

T3, T4, T4p, T5 = Q(4,5), Q(36,49), Q(85,147), Q(3,4)

@lru_cache(None)
def gb(N,d):
    K = N*delta+d+1
    return max((K+b)*(K+b+1)*
        (Q(1,4*(b+1)*(b+2))+J*tau**b)/delta**2
        for b in range(5))

@lru_cache(None)
def eb(N,b):
    low = Q(factorial(N),2**(N+1))
    for j in range(1,N+2):
        low /= b+j
    high = (rho**(N+1)-Q(2,5)**(N+1))/Q(N+1)*tau**b
    return (low+high)/delta**2

def sup(l,m,n,p):
    if m:
        return 3*T3**l*T4**(m-1)*T4p**n*T5**p
    t = l+p
    if t >= 2:
        take3, take5 = min(2,l), 2-min(2,l)
        return (Q(2)**take3*Q(3)**take5*T3**(l-take3)
                *T5**(p-take5)*T4p**n)
    if t == 1:
        return ((2*T3**(l-1)*T5**p) if l
                else (3*T5**(p-1)))*T4p**n
    return T4p**n

def ways(N,l,m,p):
    total = 0
    for t in range(l+p+1):
        multiplicity = max(0,min(l,t)-max(0,t-p)+1)
        lo = (l+p-t+m+1)//2
        hi = N-max(1,(t+m+1)//2)
        total += multiplicity*max(0,hi-lo+1)
    return total

profiles = expanded = entries = 0
Nmax = Bmax = 0
for H in range(1,68):
    for l in range(H+1):
        for m in range(H-l+1):
            for n in range(H-l-m+1):
                p = H-l-m-n
                if not p:
                    continue
                profiles += 1
                R = max(1,(l+p+1)//2+m)
                raw = Q(2)**l*Q(3)**(m+p)*T4p**n
                d = l+m+2*n+2*p
                for N in range(max(2,R),150):
                    if sup(l,m,n,p)*rho**(N-R)*gb(N,d) < 1:
                        break
                    K = N*delta+d+1
                    def good(b):
                        return raw*(K+b)*(K+b+1)*eb(N,b) < 1
                    hi = 4
                    while not good(hi):
                        hi *= 2
                    lo = 3
                    while hi > lo+1:
                        mid = (hi+lo)//2
                        if good(mid):
                            hi = mid
                        else:
                            lo = mid
                    Bcut = hi
                    w = ways(N,l,m,p)
                    if w:
                        expanded += 1
                        entries += w*Bcut
                        Nmax = max(Nmax,N)
                        Bmax = max(Bmax,Bcut)
                else:
                    raise AssertionError((H,l,m,n,p))
    if H % 5 == 0 or H == 67:
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        print(stamp,H,profiles,expanded,entries,flush=True)

assert (profiles,expanded,entries,Nmax,Bmax) == (
    916895,1738408,227336512420,71,1006)
print("PASS")

