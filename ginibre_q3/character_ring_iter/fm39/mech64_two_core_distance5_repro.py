"""FM-MECH64: exact partial closure."""
import argparse
from functools import lru_cache
from collections import defaultdict
from fractions import Fraction
from math import comb
import sympy as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--box", type=int, default=6,
                help="background exponent limit for independent checks")
args = ap.parse_args()
N,t,b,u,v = S.symbols("N t b u v")
Q = S.Rational

def choose(x,k):
    return S.prod(x-h for h in range(k))/S.factorial(k)

c = {-1:S.Integer(0), 0:S.Integer(1)}
for j in range(6):
    c[j+1] = S.expand((t*c[j]-(N-j+1)*c[j-1])/(j+1))

@lru_cache(None)
def r(j,h):
    if j < 0:
        return S.Integer(0)
    return S.expand(sum(choose(b-h,k)*c[j-2*k]
                        for k in range(j//2+1)))

def det(j,h):
    return r(j,h)**2-r(j-1,h)*r(j+1,h)

@lru_cache(None)
def pair(d,m):
    R = W = S.Integer(0)
    for h in range(d+1):
        j = d-h
        i = j-m-1
        w = 2**h*choose(b,h)
        R += w*(det(j,h)-det(i,h))
        W += w*((r(i-1,h)+r(i+1,h))*r(j,h)
                -r(i,h)*(r(j-1,h)+r(j+1,h)))
    return S.expand(R), S.expand(W)

T = {d:pair(d,d)[0] for d in range(6)}
p = {j:S.expand(c[j].subs(N,N-2*b)) for j in range(6)}
AA = {}
for d in (3,4,5):
    rem = T[d]
    A = {}
    for j in range(d,-1,-1):
        A[j] = S.factor(rem.coeff(t,2*j)*S.factorial(j)**2)
        rem = S.expand(rem-A[j]*p[j]**2)
    assert rem == 0
    AA[d] = A
    for m in range(3,d+1):
        assert S.expand(pair(d,m)[0]-T[d]
                        +(T[d-m-1] if d>m else 0)) == 0

# Exact square identities.
cert = {}
A = AA[3]
cert[3] = [A[2], A[1], A[0]-1]
for sg in (-1,1):
    rhs = ((p[3]+2*sg)**2/4 + A[2]*p[2]**2
           +A[1]*p[1]**2+A[0]-1)
    assert S.expand(rhs-pair(3,3)[0]-sg*pair(3,3)[1]) == 0

A = AA[4]
L = N+8*b-8
C2 = A[2]
C1 = A[1]-b-Q(4,5)
C0 = A[0]-1-b-L/5
cert[4] = [C2,C1,C0,A[0]-b-Q(5,4)]
for sg in (-1,1):
    rhs = (p[4]+2*sg*p[1])**2/5+L*(p[3]+2*sg)**2/20
    rhs += p[3]**2/10+b*(p[1]-sg)**2
    rhs += C2*p[2]**2+C1*p[1]**2+C0
    assert S.expand(rhs-pair(4,3)[0]-sg*pair(4,3)[1]) == 0
    rhs = (p[4]+Q(5,2)*sg)**2/5
    rhs += sum(A[j]*p[j]**2 for j in (1,2,3))
    rhs += A[0]-Q(5,4)-sg*b
    assert S.expand(rhs-pair(4,4)[0]-sg*pair(4,4)[1]) == 0

A = AA[5]
B = b*(N-4)
C2 = A[2]-Q(2,3)
C1 = A[1]-4*A[4]-Q(1,10)-B
C0 = A[0]-4*A[3]-N/2+7*b-b*N
D2 = A[2]-3*b/2
D1 = A[1]-Q(25,24)
D0 = A[0]-25*A[4]/4-Q(7,12)-b*(N-Q(3,2))
cert[5] = [
    A[4]-Q(1,10), A[3]-b, C2,C1,C0,
    A[4]-Q(1,15), D2,D1,D0,
    A[1]-b/2, A[0]-Q(3,2)-b/2
]
for sg in (-1,1):
    rhs = (p[5]+2*sg*p[2])**2/6
    rhs += (A[4]-Q(1,10))*(p[4]+2*sg*p[1])**2
    rhs += p[4]**2/10+(A[3]-b)*(p[3]+2*sg)**2+b*p[3]**2
    rhs += B*(p[1]-sg)**2+C2*p[2]**2+C1*p[1]**2+C0
    assert S.expand(rhs-pair(5,3)[0]-sg*pair(5,3)[1]) == 0

    rhs = (p[5]+Q(5,2)*sg*p[1])**2/6
    rhs += (A[4]-Q(1,15))*(p[4]+Q(5,2)*sg)**2
    rhs += p[4]**2/15+A[3]*p[3]**2+3*b*(p[2]-sg)**2/2
    rhs += D2*p[2]**2+D1*p[1]**2+D0+(1-sg)*b*(N-3)
    assert S.expand(rhs-pair(5,4)[0]-sg*pair(5,4)[1]) == 0

    rhs = (p[5]+3*sg)**2/6+b*(p[1]-sg)**2/2
    rhs += sum(A[j]*p[j]**2 for j in (2,3,4))
    rhs += (A[1]-b/2)*p[1]**2+A[0]-Q(3,2)-b/2
    assert S.expand(rhs-pair(5,5)[0]-sg*pair(5,5)[1]) == 0

# Polynomial certificates for every N >= H and b >= 1.
for d,H in ((3,4),(4,14),(5,35)):
    counts = []
    for f in cert[d]:
        z = S.Poly(S.expand(f.subs({N:H+u,b:1+v})),u,v)
        assert all(x >= 0 for x in z.coeffs())
        counts.append(len(z.terms()))
    print("large-N certificate",d,H,counts,flush=True)

# Remaining integer N,Delta: all-b tails and finite corners.
def shift(poly,B):
    return [
        sum(poly[j]*comb(j,k)*B**(j-k)
            for j in range(k,len(poly)))
        for k in range(len(poly))
    ]

for d,H,B0 in ((3,4,1),(4,14,3),(5,35,8)):
    templates = []
    for m in range(3,d+1):
        R,W = pair(d,m)
        for sg in (-1,1):
            terms = [(mon,Fraction(cc))
                     for mon,cc in S.Poly(R+sg*W,N,t,b).terms()]
            templates.append((m,sg,terms))
    tails = corners = 0
    smallmin = tailmin = None
    for NN in range(H):
        for tt in range(NN%2,NN+1,2):
            bmin = max(1,(2*d-NN+1)//2)
            for m,sg,terms in templates:
                z = [Fraction(0)]*(d+1)
                for (i,j,k),cc in terms:
                    z[k] += cc*NN**i*tt**j
                zz = shift(z,max(B0,bmin))
                assert min(zz) >= 0,(d,NN,tt,m,sg,zz)
                tails += 1
                tailmin = zz[0] if tailmin is None else min(tailmin,zz[0])
                for BB in range(bmin,B0):
                    val = sum(cc*BB**k for k,cc in enumerate(z))
                    assert val >= 0,(d,NN,tt,BB,m,sg,val)
                    corners += 1
                    smallmin = val if smallmin is None else min(smallmin,val)
    print("remaining certificates",d,"tails",tails,
          "tail minimum",tailmin,"corners",corners,
          "corner minimum",smallmin,flush=True)

# Independent SU(2) character multiplication.
def cg(i,j):
    return range(abs(i-j),i+j+1,2)

def mul(F,G):
    out = defaultdict(int)
    for (i,j),a in F.items():
        for (k,l),bb in G.items():
            for x in cg(i,k):
                for y in cg(j,l):
                    out[x,y] += a*bb
    return {key:val for key,val in out.items() if val}

S1 = {(1,0):1,(0,1):1}
D1 = {(1,0):1,(0,1):-1}
Z = {(2,0):1,(0,2):1}

def background(e,a,B):
    out = {(0,0):1}
    for fac,power in ((D1,e),(S1,a),(Z,B)):
        for _ in range(power):
            out = mul(out,fac)
    return out

def row(e,a,h):
    c = [
        sum((-1)**i*comb(e,i)*comb(a,k-i)
            for i in range(max(0,k-a),min(e,k)+1))
        for k in range(e+a+1)
    ]
    for _ in range(h):
        c = [(c[k] if k<len(c) else 0)
             +(c[k-2] if 0<=k-2<len(c) else 0)
             for k in range(len(c)+2)]
    return c

def rowpair(c,j,i):
    def at(k):
        return c[k] if 0<=k<len(c) else 0
    def DD(k):
        return at(k)**2-at(k-1)*at(k+1)
    return (DD(j)-DD(i),
            (at(i-1)+at(i+1))*at(j)
            -at(i)*(at(j-1)+at(j+1)))

def fromrows(rows,e,a,B,n,m):
    if (e+a+n+m)%2:
        return 0,0
    j = (e+a+n-m)//2
    i = (e+a+n+m)//2+1
    vv = [rowpair(rows[h],j+h,i+h) for h in range(B+1)]
    return tuple(
        sum(comb(B,h)*2**(B-h)*vv[h][q] for h in range(B+1))
        for q in (0,1)
    )

bridges = 0
for e in range(args.box+1):
    for a in range(args.box+1):
        rows = [row(e,a,h) for h in range(4)]
        for BB in range(4):
            F = background(e,a,BB)
            D = e+a+2*BB
            for n in range(D+4):
                for m in range(min(n,9)+1):
                    R = sum(F.get((j,0),0) for j in cg(n,m))
                    W = F.get((n,m),0)
                    assert (R,W) == fromrows(rows,e,a,BB,n,m)
                    if min(e,a) == 0:
                        assert R >= abs(W)
                    if (BB and m>=3 and (D+m-n)%2==0
                            and 3<=(D+m-n)//2<=5):
                        assert R >= abs(W)
                    bridges += 1
print("independent character bridges",bridges,flush=True)

# Scoped insertion obstructions.
assert rowpair(row(2,2,1),3,7) == (-1,0)
assert fromrows([row(2,2,h) for h in range(2)],
                2,2,1,3,3) == (7,0)
F = {(0,0):1,(3,3):1,(2,2):-1}
for n in range(10):
    for m in range(n+1):
        assert sum(F.get((j,0),0) for j in cg(n,m)) >= abs(F.get((n,m),0))
ZF = mul(Z,F)
assert sum(ZF.get((j,0),0) for j in cg(3,3)) == 0
assert ZF[3,3] == 2
print("scoped insertion obstructions: PASS")
print("PASS")
