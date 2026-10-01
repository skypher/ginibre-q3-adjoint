"""FM-MECH92: compatible-channel walks and the boundary normal form."""
import argparse
from collections import defaultdict
from functools import lru_cache
from math import comb

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--box", type=int, default=12)
ap.add_argument("--max-b", type=int, default=8)
ap.add_argument("--family", type=int, default=6)
args = ap.parse_args()

def C(n, j):
    return comb(n, j) if 0 <= j <= n else 0

def row(a, e):
    N = a + e
    c = [1]
    for j in range(N):
        z, rem = divmod((a-e)*c[j] -
                       (N-j+1)*(c[j-1] if j else 0), j+1)
        assert rem == 0
        c.append(z)
    return tuple(c)

def data(a, e):
    cc = row(a, e)
    N = a + e

    def c(j):
        return cc[j] if 0 <= j <= N else 0

    @lru_cache(None)
    def v(h, j):
        if h < 0:
            return 0
        if h == 0:
            return c(j)
        return v(h-1, j-1) + v(h-1, j+1)

    def D(j):
        return c(j)**2 - c(j-1)*c(j+1)

    def W(i, j):
        return c(i)*(c(j-1)+c(j+1)) - (
            c(i-1)+c(i+1))*c(j)

    @lru_cache(None)
    def T(B, j):
        if B < 0:
            return 0
        return sum(C(B,h)*2**(B-h)*(
            v(h,j)**2-v(h,j-1)*v(h,j+1))
            for h in range(B+1))

    def F(B, j):
        return T(B,j)-T(B,j+1)

    def energy(B, j):
        S = N + 2*B + 2
        H = sum(C(B,h)*2**(B-h)*(
            S*(v(h,j)**2+v(h,j-1)**2)
            -2*v(h,j-1)*((a-e)*v(h,j)
                +(B-h)*v(h+1,j)+4*h*v(h-1,j)))
            for h in range(B+1))
        return H+4*B*T(B-1,j)

    return N, c, v, D, W, T, F, energy

@lru_cache(None)
def walks(B):
    z = {(0,1): 1}
    for _ in range(B):
        nxt = defaultdict(int)
        for (i,j), w in z.items():
            nxt[i,j] += (1 if j == i+1 else 2)*w
            nxt[i-1,j-1] += w
            nxt[i-1,j+1] += w
            nxt[i+1,j+1] += w
            if j > i+1:
                nxt[i+1,j-1] += w
        z = dict(nxt)
    assert all(i < j and (j-i) % 2 and w > 0
               for (i,j),w in z.items())
    return z

def root_region(a, e, B, p):
    N = a+e
    t = min(a,e)
    x = p-2*B
    return t >= 2 and x >= 0 and (
        x*x >= 4*(t-1)*(N-t+2))

def old_band(a, e, B, p):
    S = a+e+2*B+2
    R = p*(p+2)*(S-p)*(S+p+2)
    A = abs(a-e)*(p+1)
    H = B*(S+3*p+4)
    Z = 2*R-2*A*A-H*H
    return 0 < p < S and Z >= 0 and Z*Z >= 8*A*A*H*H

def boundary(a, e, B, p, v, k):
    N = a+e
    ell = k-N
    assert ell >= 1 and p == N+2*ell
    d = B-ell
    if d < 0:
        return 0, 0
    aa = [1]
    bb = [0]
    for j in range(d+1):
        x = (a-e)*aa[j]-p*bb[j]
        y = (a-e)*bb[j]+(p+2)*aa[j]
        if j:
            x += 4*(ell+j)*aa[j-1]
            y += 4*(ell+j)*bb[j-1]
        x, rx = divmod(x, j+1)
        y, ry = divmod(y, p+j+2)
        assert rx == ry == 0
        aa.append(x)
        bb.append(y)
    sign = (-1)**e
    for j in range(d+2):
        assert aa[j] == sign*v(ell+j,k)
        assert bb[j] == sign*v(ell+j,k+1)
    value = sum(C(B,ell+j)*2**(d-j)*(
        aa[j]*bb[j+1]-aa[j+1]*bb[j]) for j in range(d+1))
    return value, d+2

# Independent Laurent-binomial construction of the walk kernel.
kernel_checks = 0
for B in range(args.max_b+1):
    z = defaultdict(int)
    for h in range(B+1):
        weight = C(B,h)*2**(B-h)
        for r in range(h+1):
            i = h-2*r
            for s in range(h+1):
                j = h-2*s+1
                w = weight*C(h,r)*C(h,s)
                if i < j:
                    z[i,j] += w
                elif j < i:
                    z[j,i] -= w
    z = {ij:w for ij,w in z.items() if w}
    assert z == walks(B)
    for ell in range(B+1):
        assert z.get((-ell,1-ell),0) >= C(B,ell)
    kernel_checks += 1
print("POSITIVE WALK KERNELS", kernel_checks, flush=True)

profiles = covered = edge = energies = channels = bbridges = 0
for a in range(args.box+1):
    for e in range(a+1):
        N,c,v,D,W,T,F,energy = data(a,e)
        for j in range((N+1)//2,N+2):
            assert D(j)-D(j+1) == W(j,j+1) >= 0
        for B in range(args.max_b+1):
            ww = walks(B)
            for k in range(N//2+1,N+B+2):
                p = 2*k-N
                val = F(B,k)
                assert val == sum(w*W(k+i,k+j)
                                  for (i,j),w in ww.items())
                profiles += 1
                h = B
                n = N-k
                assert (n+h+1)*v(h+1,k) == (
                    (a-e)*v(h,k)-p*v(h,k+1)+4*h*v(h-1,k))
                assert (k+h+2)*v(h+1,k+1) == (
                    (a-e)*v(h,k+1)+(p+2)*v(h,k)
                    +4*h*v(h-1,k+1))
                channels += 2
                loss = p*sum(C(B,h)*2**(B-h)*(
                    v(h,k-1)-v(h,k+1))**2 for h in range(B+1))
                assert energy(B,k)-energy(B,k+1) == loss
                z0 = sum(C(B,h)*2**(B-h)*v(h,k)**2
                         for h in range(B+1))
                zp = sum(C(B,h)*2**(B-h)*v(h,k+1)**2
                         for h in range(B+1))
                zm = sum(C(B,h)*2**(B-h)*v(h,k-1)**2
                         for h in range(B+1))
                assert (k+B+2)*val == (
                    (p+1)*T(B,k)+z0-zp-2*B*F(B-1,k))
                assert energy(B,k) == (
                    p*(zm-z0)+2*(k+B+1)*T(B,k)+4*B*T(B-1,k))
                energies += 1
                if k > N:
                    z, checks = boundary(a,e,B,p,v,k)
                    assert z == val
                    bbridges += checks
                if root_region(a,e,B,p):
                    assert all(W(k+i,k+j) >= 0 for i,j in ww)
                    assert val >= 0
                    covered += 1
                    if N < p <= N+2*B:
                        ell = (p-N)//2
                        assert val >= C(B,ell)
                        edge += 1
print("WALK / CONSUMER IDENTITIES", profiles,
      "CHANNEL RELATIONS", channels, flush=True)
print("ENERGY DECREMENTS", energies,
      "BOUNDARY CHANNEL BRIDGES", bbridges, flush=True)
print("ROOT REGION", covered, "NONZERO BOUNDARY", edge, flush=True)

# Definition-level Catalan moments followed by ballot conversion.
def cat(j):
    return C(2*j,j)//(j+1)

def direct(a,e,B):
    N = a+e
    degree = N+2*B
    cc = row(a,e)
    mon = [0]*(degree+1)
    for j in range(0,N+1,2):
        for l in range(B+1):
            for h in range(l+1):
                mon[N-j+2*(l-h)] += (
                    cc[j]*C(B,l)*(-2)**(B-l)*C(l,h)*cat(j//2+h))
    return [sum(mon[m]*(
        C(m,(m-p)//2)-C(m,(m-p)//2-1))
        for m in range(p,degree+1,2)) for p in range(degree+1)]

direct_checks = 0
for a in range(7):
    for e in range(7):
        N,c,v,D,W,T,F,energy = data(a,e)
        for B in range(5):
            for p,val in enumerate(direct(a,e,B)):
                if (N+p) % 2:
                    assert val == 0
                else:
                    assert val == F(B,(N+p)//2)
                direct_checks += 1
print("CATALAN BRIDGES", direct_checks, flush=True)

# An unbounded boundary family outside the earlier channel criteria.
assert -942592+519936+50128+1152 == -371376
for m in range(1,args.family+1):
    a,e,B,p = 30*m,2*m,10*m,36*m
    N,c,v,D,W,T,F,energy = data(a,e)
    k = (N+p)//2
    assert k-N == 2*m
    assert root_region(a,e,B,p)
    S=N+2*B+2
    R=p*(p+2)*(S-p)*(S+p+2)
    A=abs(a-e)*(p+1)
    H=B*(S+3*p+4)
    assert 2*R-2*A*A-H*H == (
        -942592*m**4+519936*m**3+50128*m*m+1152*m)
    assert not old_band(a,e,B,p)
    assert F(B,k) >= C(10*m,2*m) > 0
    assert boundary(a,e,B,p,v,k)[0] == F(B,k)
print("UNIFORM BOUNDARY FAMILY", args.family, "PASS", flush=True)

# An actual root-crossing control.
a=e=B=40
p=82
N,c,v,D,W,T,F,energy = data(a,e)
k=(N+p)//2
i,j=42,45
chord=W(i,j)
weight=walks(B)[i-k,j-k]
assert chord < 0 and weight > 0 and F(B,k) > 0
assert not root_region(a,e,B,p)
assert not old_band(a,e,B,p)
assert boundary(a,e,B,p,v,k)[0] == F(B,k)
print("CROSSING CONTROL", (a,e,B,p),
      "chord", chord, "weight", weight, flush=True)

# Old-band failure on (m,m,m,2m+2), m>=39.
for m in range(39,45):
    S=4*m+2
    pp=2*m+2
    R=pp*(pp+2)*(S-pp)*(S+pp+2)
    H=m*(S+3*pp+4)
    z=m-39
    assert H*H-2*R == 4*m*(z**3+81*z*z+1671*z+1239)
print("RESIDUAL FAMILY PARAMETER IDENTITY PASS", flush=True)
print("PASS", flush=True)