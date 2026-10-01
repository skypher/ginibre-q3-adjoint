import argparse
from collections import defaultdict
from fractions import Fraction as Q
from math import comb, lcm
from random import Random

ap = argparse.ArgumentParser(
    description="FM-MECH103 pair-free margin obstruction and k=2 checks")
ap.add_argument("--extended", action="store_true",
                help="also evaluate three and five copies at output 9")
args = ap.parse_args()

w = (2,2,6,5,7,4,4,1)
labs = tuple(n for n,ct in enumerate(w,1) for _ in range(ct))

def add(a,b):
    c = [0]*max(len(a),len(b))
    for i,v in enumerate(a): c[i] += v
    for i,v in enumerate(b): c[i] += v
    return c

def mul(a,b):
    c = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j] += x*y
    return c

# U_n(-2+4t).
U = [[1],[-2,4]]
for n in range(2,9):
    U.append(add(mul([-2,4],U[-1]),[-v for v in U[-2]]))

R = [1]
H = 1
for n,ct in enumerate(w,1):
    f = [-v for v in U[n]]
    f[0] += n+1
    for _ in range(ct):
        R = mul(R,f)
        H *= n+1

degree = len(R)-1
assert (sum(w),sum(labs),degree,H) == (
    31,139,139,11417363827340083200000)

K = [-20*v for v in R]
K[0] += 19*H
scale = lcm(*(comb(degree,j) for j in range(degree+1)))
bern = [
    sum(K[j]*comb(k,j)*(scale//comb(degree,j))
        for j in range(k+1))
    for k in range(degree+1)
]

def split(v):
    row = v[:]
    left = [row[0] << degree]
    right = [row[-1] << degree]
    for r in range(1,degree+1):
        row = [row[i]+row[i+1] for i in range(len(row)-1)]
        left.append(row[0] << (degree-r))
        right.append(row[-1] << (degree-r))
    return left, right[::-1]

pieces = [bern]
for _ in range(2):
    pieces = [child for parent in pieces for child in split(parent)]
assert len(pieces) == 4
assert sum(map(len,pieces)) == 560
assert all(min(v)>0 for v in pieces)
print("560 positive Bernstein coefficients: 20 R < 19 H PASS",
      flush=True)

def fusion(ns):
    st = {0:1}
    for n in ns:
        out = defaultdict(int)
        for j,v in st.items():
            for ell in range(abs(j-n),j+n+1,2):
                out[ell] += v
        st = dict(out)
    return st

# Independent two-variable character-ring multiplication.
def grow_minus(ns,p=None):
    st = {(0,0):1}
    remain = sum(ns)
    for n in ns:
        remain -= n
        out = defaultdict(int)
        for (i,j),v in st.items():
            for ell in range(abs(i-n),i+n+1,2):
                if j <= remain and (
                    p is None or abs(ell-p) <= remain
                ):
                    out[ell,j] += v
            for ell in range(abs(j-n),j+n+1,2):
                if ell <= remain and (
                    p is None or abs(i-p) <= remain
                ):
                    out[i,ell] -= v
        st = {ij:v for ij,v in out.items() if v}
    return {i:v for (i,j),v in st.items() if j == 0}

base = fusion(labs)
M = base[9]
G = grow_minus(labs,9)[9]
assert M == 24974708062978055271
assert G == 17675509457609739023
assert Q(G,M) < Q(71,100)
print("base g9/m9 =",Q(G,M),flush=True)

casimir1 = sum(n*(n+2) for n in labs)
casimir2 = sum((n*(n+2))**2 for n in labs)
assert casimir1 == 1001
assert casimir2 == 44443

for r in (1,2,3):
    st = fusion(labs*r)
    dim = H**r
    mu = 1001*r
    assert sum((j+1)*v for j,v in st.items()) == dim
    assert sum(
        (j+1)*v*j*(j+2) for j,v in st.items()
    ) == mu*dim
    second = sum(
        (j+1)*v*(j*(j+2))**2 for j,v in st.items()
    )
    assert 3*second == (5*mu*mu-2*r*casimir2)*dim
print("Casimir first/second moments PASS",flush=True)

for q in (1,2,3):
    r = q*q
    p = 40*q
    st = fusion(labs*r)
    dim = H**r
    tail = sum(
        (j+1)*v for j,v in st.items()
        if 2*j*(j+2) >= 1001*r
    )
    assert 20*tail >= 3*dim
    squares = sum(v*v for j,v in st.items() if j >= 20*q)
    assert 400*(139*r+1)**3*squares >= 9*dim*dim
    mp = fusion(labs*(2*r))[p]
    assert mp >= squares
    assert Q((p+1)**2+2012*r,(p+1)**2) < Q(903,400)

q = 40
upper = (
    Q(400,9)*(40*q+1)*(139*q*q+1)**3
    * Q(19,20)**(2*q*q)
)
assert upper < Q(1,10**50)
assert 128*Q(19,20)**114 < Q(1,2)
print("even-power lower bound and decay check PASS",flush=True)

for q in (1,2,3):
    r = q*q
    pp = 40*q+9
    oddm = fusion(labs*(2*r+1))[pp]
    evenm = fusion(labs*(2*r))[40*q]
    assert oddm >= M*evenm
    assert Q(
        (pp+1)**2+1006*(2*r+1),(pp+1)**2
    ) < Q(903,400)

q = 50
odd_upper = (
    Q(400*H,9*M)*(40*q+10)*(139*q*q+1)**3
    * Q(19,20)**(2*q*q+1)
)
assert odd_upper < Q(1,10**80)
print("full all-minus family: q=50 ratio < 10^-80 PASS",
      flush=True)

g40 = grow_minus(labs*2,40)[40]
m40 = fusion(labs*2)[40]
assert g40 > 0
print("direct g40/m40 =",Q(g40,m40),flush=True)

def at(c,j):
    return c[j] if 0 <= j < len(c) else 0

def krow(a,e):
    N = a+e
    c = [1]
    for j in range(N):
        v,r = divmod(
            (a-e)*c[-1]-(N-j+1)*(c[-2] if j else 0),j+1
        )
        assert r == 0
        c.append(v)
    return c

def td(c,j):
    return at(c,j)**2-at(c,j-1)*at(c,j+1)

def hh(c,s,ell):
    if s < ell or (s+ell)%2:
        return 0
    i,j = (s+ell)//2,(s-ell)//2
    return at(c,i)*at(c,j)-at(c,i+1)*at(c,j-1)

def k2P(a,b,n,m,d):
    c = krow(b,a+b)
    return (
        sum(td(c,d-i-j)
            for i in range(n+1) for j in range(m+1))
        - sum(hh(c,2*d-n-2*j,n) for j in range(m+1))
        - sum(hh(c,2*d-m-2*i,m) for i in range(n+1))
        + sum(hh(c,2*d-n-m,ell)
              for ell in range(abs(n-m),n+m+1,2))
    )

def coreF(c,n,m,d):
    return (
        sum(td(c,d-j)-td(c,d-n-1-j) for j in range(m+1))
        - hh(c,2*d-n,n)+hh(c,2*d-n-2*m-2,n)
        - hh(c,2*d-m,m)+hh(c,2*d-m-2*n-2,m)
        + at(c,d-m)*at(c,d-n)
        - at(c,d-m-1)*at(c,d-n-1)
        - at(c,d+1)*at(c,d-n-m-1)
        + at(c,d)*at(c,d-n-m-2)
    )

def k2F(a,b,n,m,d):
    return coreF(krow(b,a+b),n,m,d)

ns = (1,)*50+(2,)*3+(4,4)
gr = grow_minus(ns)
W = sum(ns)
co = []
acc = 0
for d in range(W//2+1):
    acc += gr.get(W-2*d,0)
    co.append(acc)

triple = co[20:23]
assert triple == [
    42580592790858697494524850,
    84746471020216909907456250,
    169457190212688866943728175
]
det = triple[1]**2-triple[0]*triple[2]
assert det == -33623261549126060189807132460564145765763518586250
assert [k2P(50,3,4,4,d) for d in (20,21,22)] == triple
assert (
    gr[22] == k2F(50,3,4,4,21)
    == triple[1]-triple[0] > 0
)
print("pair-free k=2 log-concavity failure PASS; g22",
      gr[22],flush=True)

rng = Random(10377)
checks = 0
for case in range(80):
    a = rng.randrange(10)
    b = rng.randrange(6)
    n = rng.randrange(3,11)
    m = rng.randrange(3,11)
    ns = (1,)*a+(2,)*b+(n,m)
    W = sum(ns)
    gg = grow_minus(ns)
    pp = [k2P(a,b,n,m,d) for d in range(W//2+1)]
    for d in range(W//2+1):
        f = pp[d]-(pp[d-1] if d else 0)
        assert f == gg.get(W-2*d,0) == k2F(a,b,n,m,d)
        p = W-2*d
        psign = (-1)**len(ns)
        if p >= max(n,m) and not (
            psign == 1 and p in (n,m)
        ):
            assert f >= 0
        checks += 1
assert checks == 966
print("80 k=2 profiles; 966 coefficient bridges PASS",flush=True)

# Free recurrence states need not give a PSD quadratic form.
c = [Q(0)]*67
c[7],c[8] = Q(1),Q(-2)
for j in range(7,0,-1):
    c[j-1] = (-64*c[j]-(j+1)*c[j+1])/(67-j)
c[9] = (-64*c[8]-59*c[7])/9
assert coreF(c,3,3,8) == -Q(
    8681792848128739,586315037107800)

actual = krow(1,65)
assert Q(actual[8],actual[7]) == -Q(1475,208)
g56 = grow_minus((1,)*64+(2,3,3),56)[56]
assert g56 == coreF(actual,3,3,8) > 0
print("free-state Gram failure PASS; actual g56",
      g56,flush=True)

if args.extended:
    for copies,bound in ((3,Q(1,2)),(5,Q(1,3))):
        gg = grow_minus(labs*copies,9)[9]
        mm = fusion(labs*copies)[9]
        assert 0 < Q(gg,mm) < bound
        print("copies",copies,"ratio",Q(gg,mm),flush=True)

print("ALL CHECKS PASS",flush=True)