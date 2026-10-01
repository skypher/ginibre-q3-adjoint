"""FM-MECH60 verifier. Run from the repository root with python3 -u.
The C<4096 box is part of the proof; the --audit-max-n census is vetting.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import sympy as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--audit-max-n", type=int, default=120)
args = ap.parse_args()

# Local metric, fixed energy, and Sonin remainder.
sig, X, de, om, p, q = S.symbols("sig X de om p q")
C = sig**2-de**2
nxt = (2*de*p-(sig-X)*q)/(sig+X)
nxt2 = (2*de*nxt-(sig-X-2)*p)/(sig+X+2)
B = q+nxt
D = p*p-q*nxt
Dn = nxt*nxt-p*nxt2
H = sig*(p*p+q*q)-2*de*p*q
Hn = sig*(nxt*nxt+p*p)-2*de*nxt*p
Q = S.Matrix([
    [1+de**2/X**2, -de*sig/(2*X**2)],
    [-de*sig/(2*X**2), -S.Rational(1,4)+sig**2/(4*X**2)]
])
assert S.cancel(D-(S.Matrix([p,B]).T*Q*S.Matrix([p,B]))[0]) == 0
assert S.factor(Q.det()-(C-X**2)/(4*X**2)) == 0
assert S.cancel(H-Hn-4*X*(sig*q-de*p)**2/(sig+X)**2) == 0
assert S.cancel((sig+X)*D-H-X*(p*p-q*q)) == 0
assert S.expand(H*H-C*(p*p-q*q)**2
                -(de*(p*p+q*q)-2*sig*p*q)**2) == 0
rate = (sig-X)*(om+X+2)/((sig+X+2)*(om+X))
rem = 2*((sig-om)*nxt*nxt-2*de*nxt*p+(sig+om)*p*p)/(
    (sig+X+2)*(om+X))
assert S.cancel(rate*D-Dn-rem) == 0

# Modified energy; coordinates here are scaled by omega.
V = C+2*sig-X*(X+2)
u0 = sig*p-de*B/2
v0 = (sig*B/2-de*p)/X
u1 = sig*nxt-de*(p+nxt2)/2
v1 = (sig*(p+nxt2)/2-de*nxt)/(X+2)
assert S.cancel(u1-(de*u0-V*v0)/(sig+X+2)) == 0
assert S.cancel(v1-(u0+de*v0)/(sig+X+2)) == 0
assert S.cancel(u0*u0+(C-X*X)*v0*v0-C*D) == 0
E0 = u0*u0+V*v0*v0
E1 = u1*u1+(V-4*(X+2))*v1*v1
assert S.cancel((sig-X)*E0/(sig+X+2)-E1-4*(X+2)*v1*v1) == 0
print("METRIC, SONIN, BINOMIAL ENERGY IDENTITIES PASS")

# Exact lower Taylor bounds for the exponentials used in the proof.
assert sum(F(45,32)**h/factorial(h) for h in range(4)) > F(18,5)
assert sum(F(15,4)**h/factorial(h) for h in range(13)) > 42
assert F(4,3)+F(32,7) < F(15,2)
assert F(29,30)**2-F(1089,42)*F(1,30) == F(221,3150)
assert sum(F(1107,640)**h/factorial(h) for h in range(5)) > F(845,168)
assert sum(F(2907,640)**h/factorial(h) for h in range(16)) > 90
for y in (F(5,8), F(29,32)):
    assert 15*y**3-16*y+3 < 0
for y in (F(7,8), F(29,32)):
    assert 15*y**3-15*y+2 < 0
assert F(331,336)**2-4*F(33,32)**2*F(25,126) == F(28547,225792)
print("EXPONENTIAL AND CHECKPOINT CONSTANTS PASS")

def data(a, e):
    N = a+e
    si, dd = N+2, a-e
    cc = si*si-dd*dd
    c = [1]
    for k in range(N):
        z, r = divmod(dd*c[k]-(N-k+1)*(c[k-1] if k else 0), k+1)
        assert r == 0
        c.append(z)
    v = lambda k: c[k] if 0 <= k <= N else 0
    b = [v(k-1)+v(k+1) for k in range(N+2)]
    d = [v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    h = [si*(v(k)**2+v(k-1)**2)-2*dd*v(k)*v(k-1)
         for k in range(N+2)]
    return N, si, dd, cc, v, b, d, h

def corridor(x, y, c):
    return (16*x*x <= c and 4*y*y >= c or
            25*x*x <= 4*c and 64*y*y >= 25*c)

# Fixed finite remainder completing both uniform regions.
rows = pairs = 0
least = None
for e in range(3,31):
    for a in range(e+2,256):
        si = a+e+2
        if 4*(a-e) <= si or 4*(a+1)*(e+1) >= 4096:
            continue
        N, si, dd, cc, v, b, d, h = data(a,e)
        rows += 1
        for j in range(N//2+1,N-3):
            xx = 2*j-N
            if 25*xx*xx > 4*cc:
                continue
            for i in range(j+4,N+2):
                yy = 2*i-N
                if not corridor(xx,yy,cc):
                    continue
                w = v(j)*b[i]-b[j]*v(i)
                slack = d[j]-d[i]-abs(w)
                assert slack >= 0, (a,e,j,i)
                pairs += 1
                item = slack, (a,e,j,i)
                if least is None or item < least:
                    least = item
assert (rows,pairs,least) == (1494,490037,(67,(6,3,5,9)))
print("FINITE CORRIDOR", rows, pairs, "MINIMUM", least)

# Independent definition-level Catalan reduction.
@lru_cache(None)
def moment(power, label):
    if (power+label) % 2:
        return 0
    return sum(
        (-1)**h*comb(label-h,h)*
        comb(power+label-2*h,(power+label-2*h)//2)//
        ((power+label-2*h)//2+1)
        for h in range(label//2+1)
    )

@lru_cache(None)
def kernel(e, a, p, q):
    return sum(
        (-1)**h*comb(e,h)*comb(a,k)*
        moment(e+a-h-k,p)*moment(h+k,q)
        for h in range(e+1) for k in range(a+1)
    )

bridges = 0
for N in range(8,23):
    for e in range(3,(N-2)//2+1):
        a = N-e
        if 4*(a-e) <= N+2:
            continue
        N, si, dd, cc, v, b, d, h = data(a,e)
        for j in range(N//2+1,N-3):
            if 25*(2*j-N)**2 > 4*cc:
                continue
            for i in range(j+4,N+2):
                if not corridor(2*j-N,2*i-N,cc):
                    continue
                p, q = i+j-N-1, i-j-1
                assert kernel(e,a,p,q) == v(j)*b[i]-b[j]*v(i)
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2)) == d[j]-d[i]
                bridges += 1
print("CATALAN BRIDGES", bridges)

def nonnegative_radical(A, B, cc):
    """A+B*sqrt(cc)>=0, using exact arithmetic."""
    if B == 0:
        return A >= 0
    if B > 0:
        return A >= 0 or B*B*cc >= A*A
    return A >= 0 and A*A >= B*B*cc

# Diagnostic census; this is not an unbounded proof obligation.
total = pcert = bcert = actual = 0
misses = []
for N in range(8,args.audit_max_n+1):
    for e in range(3,(N-2)//2+1):
        a = N-e
        if 4*(a-e) <= N+2:
            continue
        N, si, dd, cc, v, b, d, h = data(a,e)
        binom = [comb(si,k+1) for k in range(N+1)]
        for j in range(N//2+1,N-3):
            xx = 2*j-N
            if 2*xx > si or xx*xx >= cc:
                continue
            crossed = False
            for i in range(j+1,N+1):
                yy = 2*i-N
                if yy*yy >= cc:
                    break
                w = v(j)*b[i]-b[j]*v(i)
                if w < 0 or (w == 0 and v(j)*v(i)+b[j]*b[i] < 0):
                    crossed = True
                if not crossed or i-j < 4:
                    continue
                total += 1

                A = (cc-xx*yy)*binom[j]-(cc+yy*yy)*binom[i]
                BB = (xx-yy)*binom[j]-2*yy*binom[i]
                ok = nonnegative_radical(A,BB,cc)
                pcert += ok

                P = cc-xx*xx
                V = P+2*(si-xx)
                okb = nonnegative_radical(
                    -yy*(P*binom[j]+V*binom[i]),
                    P*binom[j]-V*binom[i], cc)
                bcert += okb

                local = ((cc-yy*yy)*(d[j]-d[i])**2
                         >= 4*yy*yy*d[j]*d[i])
                actual += local
                assert not (ok or okb) or local
                assert d[j]-d[i] >= abs(w)
                if not ok:
                    misses.append((a,e,j,i))
if args.audit_max_n == 120:
    assert (total,pcert,bcert,actual) == (
        998772,998765,998772,998772)
print("INNER LONG AUDIT", total, "SONIN", pcert,
      "BINOMIAL", bcert, "ACTUAL METRIC", actual)
print("PARAMETER MISSES", misses)

# Earlier Sonin-envelope failure, repaired by the new energy.
N, si, dd, cc, v, b, d, h = data(37,4)
j, i = 21, 25
xx, yy = 2*j-N, 2*i-N
z = F(comb(si,i+1),comb(si,j+1))
A = cc-xx*yy-(cc+yy*yy)*z
BB = xx-yy-2*yy*z
assert z == F(1197,2990)
assert A == F(1238813,2990) and BB == -F(22733,1495)
assert A > 0 and BB < 0 and A*A < BB*BB*cc
assert v(j)*b[i]-b[j]*v(i) < 0
assert (cc-yy*yy)*(d[j]-d[i])**2 >= 4*yy*yy*d[j]*d[i]
P = cc-xx*xx
V = P+2*(si-xx)
assert nonnegative_radical(-yy*(P+V*z),P-V*z,cc)
print("SONIN FAILURE REPAIRED", A, BB, "sqrt", cc,
      "ACTUAL E SLACK", d[j]-d[i]-abs(v(j)*b[i]-b[j]*v(i)))

# Near-turning failure of the local metric, covered by C2.
N, si, dd, cc, v, b, d, h = data(4000,3)
j, i = 2034, 2128
yy = 2*i-N
drop = d[j]-d[i]
local = F((cc-yy*yy)*drop*drop,4*yy*yy*d[j]*d[i])
assert local < 1
assert cc*(si+yy)**2*drop*drop >= 4*yy*yy*h[j]*h[i]
assert drop >= abs(v(j)*b[i]-b[j]*v(i))
assert corridor(2*j-N,yy,cc) and 16*(2*j-N)**2 > cc
print("TURNING METRIC FAILURE",
      (local.numerator*10**6//local.denominator,
       local.numerator*10**6//local.denominator+1),
      "OVER 10^6; CORRIDOR, J, E PASS")
print("PASS")

