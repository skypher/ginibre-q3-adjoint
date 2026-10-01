import argparse
from math import comb
from functools import lru_cache
from fractions import Fraction as Q
import sympy as S

argparse.ArgumentParser(description="FM-MECH65 exact cone tests").parse_args()

def C(n,k): return comb(n,k) if 0<=k<=n else 0
def cat(k): return C(2*k,k)//(k+1)
def trim(v):
    v=list(v)
    while v and v[-1]==0: v.pop()
    return tuple(v)
def at(v,p): return v[p] if 0<=p<len(v) else 0

@lru_cache(None)
def moment(e,a,b,k=0):
    T=e+a; D=T+2*b
    row=[sum((-1)**h*C(e,h)*C(a,j-h)
             for h in range(e+1)) for j in range(T+1)]
    mon=[0]*(D+1)
    for j in range(k%2,T+1,2):
        for u in range(b+1):
            for v in range(b-u+1):
                mon[T-j+2*u]+=(
                    row[j]*C(b,u)*C(b-u,v)*(-2)**(b-u-v)
                    *cat((j+k)//2+v))
    return trim([
        sum(mon[n]*(C(n,(n-p)//2)-C(n,(n-p)//2-1))
            for n in range(p,D+1,2))
        for p in range(D+1)])

def X(v):
    out=[0]*(len(v)+1)
    for i,c in enumerate(v):
        out[i+1]+=c
        if i: out[i-1]+=c
    return trim(out)
def add(*vs):
    n=max((len(v) for _,v in vs),default=0)
    return trim([sum(k*at(v,i) for k,v in vs) for i in range(n)])
def transfer(u,v,nextu):
    return (
        add((1,X(u)),(-1,v)),
        add((1,X(v)),(1,X(X(u))),(-2,u),(-1,nextu)))

cases=[
    (e,T-e,b,T+2*b-2*d)
    for T in (13,14) for b in range(3,10) for d in (6,7)
    if T+2*b-2*d>=7 for e in range(2,T-1)]
cases += [(2,2,15,22),(23,2,3,19)]
for e,a,b,p in cases:
    u,v=moment(e-1,a,b),moment(e-1,a,b,1)
    uu,vv=transfer(u,v,moment(e-1,a,b+1))
    assert uu==moment(e,a,b) and vv==moment(e,a,b,1)
    assert min(uu,default=0)>=0
    A,B=at(X(u),p),at(v,p)
    assert A-B==at(uu,p)
    assert A+B==at(moment(e-1,a+1,b),p)
    assert A>=abs(B)
assert len(cases)==275
u,v=moment(5,8,3),moment(5,8,3,1)
assert (at(X(u),8),at(v,8))==(3136,1016)
print("275 full-vector transfer checks passed; target =",3136-1016)

for n in range(2,21):
    vals=[at(moment(0,n,b),n) for b in range(3)]
    assert vals==[
        1,n*(n+1)//2,(n**4+17*n*n-6*n+24)//12]
    t=n-2
    assert vals[0]*vals[2]-vals[1]**2 == -Q(
        t**4+11*t**3+35*t*t+43*t+6,6)

vals=[at(moment(6,8,b),8) for b in range(3,8)]
H=S.Matrix(3,3,lambda i,j:vals[i+j])
q=S.Matrix([17,-12,2])
assert vals==[2120,5205,13255,35370,99330]
assert H.det()==-340496125
assert (q.T*H*q)[0]==-1340

blocks=[]
for b in range(3,8):
    u,v=moment(5,8,b),moment(5,8,b,1)
    A,B=at(X(u),8),at(v,8)
    blocks.append(S.Matrix([[A,B],[B,A]]))
HH=S.Matrix(6,6,lambda i,j:blocks[i//2+j//2][i%2,j%2])
qq=S.Matrix([z for x in [17,-12,2] for z in [x,-x]])
assert (qq.T*HH*qq)[0]==-2680
print("Residual Hankel determinant:",H.det())
print("Coupled Hankel negative quadratic:",-2680)

K0=S.Matrix([[at(moment(0,2,b),p) for p in (0,2)]
             for b in (0,1)])
K1=S.Matrix([[at(moment(0,2,b),p) for p in (0,2)]
             for b in (1,2)])
assert K0==S.Matrix([[2,1],[2,3]]) and K0.det()==4
assert K1==S.Matrix([[2,3],[6,8]]) and K1.det()==-2

# Exact fixed-b polynomials and their leading coefficients.
n=S.symbols("n")
def binpoly(z,k):
    return (S.prod(z-i for i in range(k))/S.factorial(k)
            if k>=0 else S.Integer(0))
def cpoly(b):
    ans=0
    for j in range(b+1):
        for u in range(j,b+1):
            h=u-j
            for v in range(b-u+1):
                ans+=(
                    binpoly(n,2*j)*C(b,u)*C(b-u,v)
                    *(-2)**(b-u-v)*cat(j+v)
                    *(binpoly(n+2*h,h)-binpoly(n+2*h,h-1)))
    return S.Poly(S.expand(ans),n)
ps=[cpoly(b) for b in range(6)]
for b,P in enumerate(ps):
    assert P.degree()==2*b
    assert P.LC()==1/(S.factorial(b)*S.factorial(b+1))
    for a in range(2,9):
        assert P.eval(a)==at(moment(0,a,b),a)
for s in range(4):
    D=S.Poly(ps[s].as_expr()*ps[s+2].as_expr()
             -ps[s+1].as_expr()**2,n)
    assert D.degree()==4*s+4
    assert D.LC()==-S.Rational(
        2,(s+3)*(S.factorial(s+1)*S.factorial(s+2))**2)
assert [at(moment(0,7,b),7) for b in (3,4,5)]==[
    1808,10904,62528]
assert 1808*62528-10904**2==-5846592

# Failure of the two-component order cone.
u,v=(1,),()
tu,tv=transfer(u,v,u)
assert tu==(0,1) and tv==(-2,0,1)
assert (at(X(tu),0),at(tv,0))==(1,-2)
assert S.Matrix([[1,-2],[-2,1]]).det()==-3
uu,vv=transfer(tu,tv,tu)
assert uu==(3,)
rr,_=transfer(uu,vv,uu)
assert rr==(0,5,0,-2)

# Signed-character cone: exact generator preservation.
def bimul(f,n,axis,sg=1):
    out={}
    for (p,q),v in f.items():
        r=p if axis==0 else q
        for t in range(abs(r-n),r+n+1,2):
            key=(t,q) if axis==0 else (p,t)
            out[key]=out.get(key,0)+sg*v
    return out
def bisum(f,g):
    out=dict(f)
    for k,v in g.items(): out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}
for p in range(7):
    for q in range(7):
        f={(p,q):(-1)**q}
        D=bisum(bimul(f,1,0),bimul(f,1,1,-1))
        Z=bisum(bimul(f,2,0),bimul(f,2,1))
        assert all((-1)**j*v>=0 for (i,j),v in D.items())
        assert all((-1)**j*v>=0 for (i,j),v in Z.items())
seed2={(2,0):1,(0,2):1,(1,1):2,(0,0):2}
assert -seed2[1,1]==-2
print("Uniform seed formulas, tail formulas, TP/order obstructions,")
print("and signed-character generator checks: PASS")
