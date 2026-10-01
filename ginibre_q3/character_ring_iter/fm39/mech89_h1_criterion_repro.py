"""FM-MECH89: channel constraints and a uniform central region."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from random import Random
import sympy as sp

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--cases", type=int, default=400)
ap.add_argument("--family", type=int, default=12)
args = ap.parse_args()

def choose(n, j):
    return comb(n, j) if 0 <= j <= n else 0

@lru_cache(None)
def row(a, e):
    N = a+e
    c = [1]
    for k in range(N):
        v, rem = divmod(
            (a-e)*c[k]-(N-k+1)*(c[k-1] if k else 0), k+1)
        assert rem == 0
        c.append(v)
    return tuple(c)

def data(a, e, B, p):
    N = a+e
    assert (N+p) % 2 == 0
    k = (N+p)//2
    cc = row(a, e)

    def c(l):
        return cc[l] if 0 <= l <= N else 0

    @lru_cache(None)
    def v(h, l):
        if h < 0:
            return 0
        return sum(comb(h, i)*c(l+h-2*i) for i in range(h+1))

    def T(b, l):
        if b < 0:
            return 0
        return sum(
            comb(b, h)*2**(b-h)*
            (v(h, l)**2-v(h, l-1)*v(h, l+1))
            for h in range(b+1))

    def F(b, l):
        return T(b, l)-T(b, l+1)

    def norm(l):
        return sum(comb(B, h)*2**(B-h)*v(h, l)**2
                   for h in range(B+1))

    def dotA(l, m):
        return sum(
            comb(B, h)*2**(B-h)*v(h, l)*
            ((a-e)*v(h, m)+(B-h)*v(h+1, m)+4*h*v(h-1, m))
            for h in range(B+1))

    return N, k, v, T, F, norm, dotA

def criterion(n, p, B, delta):
    if n < 2*B:
        return False
    q = p*(p+2)
    L = n+1-2*B
    A = abs(delta)*(p+1)
    C2 = 32*B*B*(q+delta*delta)
    Z = 4*q*L*L-A*A-C2
    return Z >= 0 and Z*Z >= 4*A*A*C2

def old_band(N, p, B, delta):
    S = N+2*B+2
    R = p*(p+2)*(S-p)*(S+p+2)
    A = abs(delta)*(p+1)
    C = B*(S+3*p+4)
    Z = 2*R-2*A*A-C*C
    return 0 < p < S and Z >= 0 and Z*Z >= 8*A*A*C*C

def gram(n, p, B, delta, v, k):
    lam = delta*delta+p*(p+2)
    out = Q(0)
    for j in range(B//2+1):
        w = comb(B, 2*j)*2**(B-2*j)
        o = (comb(B, 2*j+1)*2**(B-2*j-1)
             if 2*j < B else 0)
        a = n+2*j+1
        b = n+p+2*j+2
        C = Q(w)-Q(4*(2*j+1)*o, a+1)
        D = Q(w)-Q(4*(2*j+1)*o, b+1)
        r = v(2*j+1, k+1)
        s = -v(2*j+1, k)
        rm = v(2*j-1, k+1)
        sm = -v(2*j-1, k)
        out += (
            C*p*b*r*r+D*(p+2)*a*s*s
            +delta*(D*b-C*a)*r*s
            -8*j*C*p*r*rm-8*j*D*(p+2)*s*sm
            +8*j*delta*(C*r*sm-D*s*rm)
            +lam*o*(
                Q(p, a+1)*r*r+Q(p+2, b+1)*s*s
                +Q(delta*(p+1), (a+1)*(b+1))*r*s))
    return out

# Universal single-pair identity, with denominators cleared.
n, p, j, d, w, o, r, s, rm, sm = sp.symbols(
    "n p j d w o r s rm sm")
aa = n+2*j+1
bb = n+p+2*j+2
lam = d*d+p*(p+2)
xnum = -d*(aa*s-8*j*sm)+p*(bb*r-8*j*rm)
ynum = (p+2)*(aa*s-8*j*sm)+d*(bb*r-8*j*rm)
xnextnum = (-d*s-p*r)*lam+4*(2*j+1)*xnum
ynextnum = (-(p+2)*s+d*r)*lam+4*(2*j+1)*ynum
Cnum = w*(aa+1)-4*(2*j+1)*o
Dnum = w*(bb+1)-4*(2*j+1)*o
left = (
    w*(aa+1)*(bb+1)*(xnum*r+s*ynum)
    -o*((aa+1)*s*ynextnum+(bb+1)*xnextnum*r))
right = (
    Cnum*(bb+1)*p*bb*r*r
    +Dnum*(aa+1)*(p+2)*aa*s*s
    +d*(Dnum*(aa+1)*bb-Cnum*(bb+1)*aa)*r*s
    -8*j*Cnum*(bb+1)*p*r*rm
    -8*j*Dnum*(aa+1)*(p+2)*s*sm
    +8*j*d*(Cnum*(bb+1)*r*sm-Dnum*(aa+1)*s*rm)
    +lam*o*(
        p*(bb+1)*r*r+(p+2)*(aa+1)*s*s
        +d*(p+1)*r*s))
assert sp.expand(left-right) == 0

rho = sp.symbols("rho")
assert sp.expand(
    bb*((bb+1)-rho)*(aa+1)
    -aa*((aa+1)-rho)*(bb+1)
    -(p+1)*((aa+1)*(bb+1)-rho)) == 0
print("UNIVERSAL BLOCK-GRAM ALGEBRA PASS", flush=True)

rng = Random(8901)
profiles = []
for _ in range(args.cases):
    B = rng.randrange(1, 25)
    n0 = rng.randrange(0, 10*B+1)
    p0 = rng.randrange(1, 41)
    N = 2*n0+p0
    delta = rng.randrange(-N, N+1, 2)
    profiles.append(((N+delta)//2, (N-delta)//2, B, p0))

# Both proved uniform subregions.
for _ in range(80):
    B = rng.randrange(1, 33)
    p0 = 2*rng.randrange(1, 25)
    n0 = 5*B+rng.randrange(0, 4*B+1)
    N = 2*n0+p0
    profiles.append((N//2, N//2, B, p0))

for _ in range(80):
    B = rng.randrange(1, 33)
    p0 = rng.randrange(1, 49)
    n0 = 8*B+rng.randrange(0, 4*B+1)
    N = 2*n0+p0
    limit = min(p0, n0//4)
    delta = rng.randrange(-limit, limit+1)
    if (N+delta) % 2:
        delta = delta+1 if delta < limit else delta-1
    profiles.append(((N+delta)//2, (N-delta)//2, B, p0))

identities = channels = covered = balanced = near = outside_old = 0
for a, e, B, p0 in profiles:
    N, k, v, T, F, norm, dotA = data(a, e, B, p0)
    n0 = N-k
    delta = a-e
    value = F(B, k)

    for h in range(B+2):
        assert (n0+h+1)*v(h+1, k) == (
            delta*v(h, k)-p0*v(h, k+1)+4*h*v(h-1, k))
        assert (k+h+2)*v(h+1, k+1) == (
            delta*v(h, k+1)+(p0+2)*v(h, k)
            +4*h*v(h-1, k+1))
        channels += 1

    assert gram(n0, p0, B, delta, v, k) == (
        delta*delta+p0*(p0+2))*value

    S = N+2*B+2
    t = k+B+1
    u = N+B-k+1
    left = dotA(k, k+1)
    right = dotA(k+1, k)
    assert left-right == -2*B*F(B-1, k)
    assert 2*u*(t+1)*value == (
        2*u*(p0+2)*norm(k)+2*p0*(t+1)*norm(k+1)
        -(p0+1)*(left+right)-2*B*(S+1)*F(B-1, k))
    identities += 1

    if delta == 0 and n0 >= 5*B:
        assert criterion(n0, p0, B, delta)
        balanced += 1
    if n0 >= 8*B and abs(delta) <= min(p0, n0//4):
        assert criterion(n0, p0, B, delta)
        near += 1
    if criterion(n0, p0, B, delta):
        assert value >= 0
        covered += 1
        if not old_band(N, p0, B, delta):
            outside_old += 1

print("ACTUAL PROFILES", identities,
      "CHANNEL RELATIONS", channels, flush=True)
print("CRITERION", covered, "BALANCED", balanced,
      "NEAR-BALANCED", near,
      "OUTSIDE MECH86 BAND", outside_old, flush=True)

# Independent Catalan moments and ballot conversion.
def cat(h):
    return comb(2*h, h)//(h+1)

def direct(a, e, B):
    N = a+e
    degree = N+2*B
    mon = [0]*(degree+1)
    cc = row(a, e)
    for j0 in range(0, N+1, 2):
        for l in range(B+1):
            for h in range(l+1):
                mon[N-j0+2*(l-h)] += (
                    cc[j0]*comb(B, l)*(-2)**(B-l)
                    *comb(l, h)*cat(j0//2+h))
    return [
        sum(mon[d0]*(
            choose(d0, (d0-p0)//2)
            -choose(d0, (d0-p0)//2-1))
            for d0 in range(p0, degree+1, 2))
        for p0 in range(degree+1)]

bridges = 0
for a in range(7):
    for e in range(7):
        for B in range(5):
            vals = direct(a, e, B)
            for p0, value in enumerate(vals):
                if (a+e+p0) % 2:
                    assert value == 0
                else:
                    N, k, v, T, F, norm, dotA = data(a, e, B, p0)
                    assert F(B, k) == value
                bridges += 1
print("CATALAN BRIDGES", bridges, flush=True)

# The uniform family lies outside the MECH86 band identically.
m = sp.symbols("m", positive=True, integer=True)
NN, BB, pp = 84*m, 8*m, 4*m
SS = NN+2*BB+2
old_difference = sp.expand(
    BB**2*(SS+3*pp+4)**2
    -2*pp*(pp+2)*(SS-pp)*(SS+pp+2))
assert all(c > 0 for c in sp.Poly(
    old_difference.subs(m, m+1), m).coeffs())

for m0 in range(1, args.family+1):
    a = e = 42*m0
    B = 8*m0
    p0 = 4*m0
    N, k, v, T, F, norm, dotA = data(a, e, B, p0)
    assert N-k == 5*B
    assert criterion(N-k, p0, B, 0)
    assert not old_band(N, p0, B, 0)
    assert F(B, k) > 0
print("UNIFORM FAMILY", args.family,
      "PASS; outside old band identically", flush=True)

# The block form is not PSD on independent odd-channel data.
n0, p0, B = 108, 12, 108
a0 = n0+1
b0 = n0+p0+2
C0 = Q(1)-Q(2*B, a0+1)
free_diag = b0*C0+Q(B*p0*(p0+2), 2*(a0+1))
assert free_diag == Q(-386, 11)

x0 = Q(122, 14)
x1_forward = (4*x0-12)/110
x1_next_pair = Q(-8, 14)
assert x1_forward == Q(16, 77)
assert x1_next_pair == Q(-4, 7)
assert x1_forward != x1_next_pair

N, k, v, T, F, norm, dotA = data(114, 114, 108, 12)
actual = F(108, k)
assert actual > 0
print("FREE-CHANNEL GRAM CONTROL", free_diag,
      "times 2^108", flush=True)
print("MISSING COMPATIBILITY CONTROL",
      x1_forward, x1_next_pair, flush=True)
print("ACTUAL CONTROL VALUE", actual, flush=True)
print("PASS", flush=True)