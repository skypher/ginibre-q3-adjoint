import argparse
argparse.ArgumentParser(
    description="FM-MECH88 exact verifier; memory only"
).parse_args()
from math import comb
from functools import lru_cache
from random import Random
import sympy as S

def mul(a,b,cap):
    out={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            if i+k+2*(j+l)<=cap:
                ij=(i+k,j+l)
                out[ij]=out.get(ij,0)+v*w
    return {m:v for m,v in out.items() if v}

def power(a,n,cap):
    out={(0,0):1}
    while n:
        if n&1:
            out=mul(out,a,cap)
        n//=2
        if n:
            a=mul(a,a,cap)
    return out

def block(n,e,cap):
    f={(0,j):1 for j in range(min(n,cap//2)+1)}
    if n<=cap:
        for j in range(n//2+1):
            ij=(n-2*j,j)
            f[ij]=f.get(ij,0)+e*(-1)**j*comb(n-j,j)
    return {ij:v for ij,v in f.items() if v}

def poly(cs,d):
    out={(0,0):1}
    for n,e,c in cs:
        out=mul(out,power(block(n,e,2*d),c,2*d),2*d)
    return out

def gram(cs,d):
    pp=poly(cs,d)
    def entry(i,j):
        out=0
        for h in range(min(i,j)+1):
            r=i+j-2*h
            out+=comb(r,i-h)*pp.get((r,h),0)*S.Rational(2,r+2)
        return out
    return S.Matrix(d+1,d+1,entry)

def pvals(cs,d):
    pp=poly(cs,d)
    return [
        sum(v*(comb(i,i//2)//(i//2+1))
            for (i,j),v in pp.items()
            if i%2==0 and i//2+j==r)
        for r in range(d+1)
    ]

def difference(G):
    return S.Matrix(
        G.rows,G.cols,
        lambda i,j:G[i,j]-(G[i-1,j-1] if min(i,j) else 0)
    )

def psd(A):
    A=S.MutableDenseMatrix(A)
    for j in range(A.rows):
        p=A[j,j]
        if p<0:
            return False,j,p
        if p==0:
            if any(A[j,i] for i in range(j+1,A.rows)):
                return False,j,'zero-cross'
            continue
        for i in range(j+1,A.rows):
            for k in range(i,A.rows):
                v=A[i,k]-A[i,j]*A[j,k]/p
                A[i,k]=A[k,i]=v
    return True

# Independent Clebsch-Gordan coefficient calculation.
def build(cs,cap):
    out={(0,0):1}
    for n,e,c in cs:
        for _ in range(c):
            new={}
            for (a,j),v in out.items():
                for r in range(n+1):
                    aa=a+2*r
                    if aa<=cap:
                        new[aa,j]=new.get((aa,j),0)+v
                if a+n<=cap:
                    for ell in range(abs(j-n),j+n+1,2):
                        new[a+n,ell]=new.get((a+n,ell),0)+e*v
            out={ij:v for ij,v in new.items() if v}
    return out

def F(cs,d):
    T=build(cs,2*d)
    return T.get((2*d,0),0)-T.get((2*d-2,0),0)

def neutral_formula(bg,n,l):
    d=n+l
    T=build(bg,2*d)
    P=[T.get((2*j,0),0) for j in range(d+1)]
    def tor(j):
        return sum(v for (a,k),v in T.items()
                   if a==2*j and k%2==0)
    return sum(P)-2*sum(P[:l])-tor(l)+(tor(l-1) if l else 0)

rng=Random(8806)
done=0
for case in range(100):
    bg=[(rng.randrange(1,8),rng.choice((-1,1)),rng.randrange(1,4))
        for _ in range(rng.randrange(2,5))]
    D=sum(n*c for n,e,c in bg)
    M=max(n for n,e,c in bg)
    l=rng.randrange(1,max(2,min(5,(D-M)//2+1)))
    p=D-2*l
    if p<max(M,l):
        continue
    n=rng.randrange(l,min(p,l+8)+1)
    d=n+l
    old=F(bg+[(n,1,1),(n,-1,1)],d)
    new=F(bg+[(l,1,1),(l,-1,1)],2*l)
    pp=pvals(bg,d)
    assert old==neutral_formula(bg,n,l)
    assert old-new==sum(pp[2*l+1:d+1])>=n-l
    done+=1
assert done==100
print("neutral reflection: 100 admissible profiles PASS",flush=True)

for n in range(3,15):
    bg=[(1,1,n),(2,-1,1)]
    assert F(bg+[(n,1,1),(n,-1,1)],n)==sum(pvals(bg,n)[1:])
print("zero-deficit boundary PASS",flush=True)

def endpoint(cs,d):
    row=[1]+[0]*d
    for n,e,c in cs:
        for _ in range(c):
            for j in range(d,n-1,-1):
                row[j]+=e*row[j-n]
    return row[d]

def endpoint_cert(d,A,e):
    return e*A>=0 or abs(A)<=1 or abs(A)>=d

rng=Random(8808)
checked=certified=0
for case in range(100):
    d=rng.randrange(3,14)
    bg=[(rng.randrange(1,d+1),rng.choice((-1,1)),rng.randrange(1,4))
        for _ in range(rng.randrange(2,6))]
    D=sum(n*c for n,e,c in bg)
    if D<2*d:
        continue
    T=build(bg,2*d)
    p=T.get((2*d,0),0)
    A=endpoint(bg,d)
    assert (d+1)*p>=A*A and p>=1
    for e in (-1,1):
        f=F(bg+[(d,e,1)],d)
        assert f==p+e*A
        if endpoint_cert(d,A,e):
            assert f>=0
            certified+=1
    checked+=1
print("endpoint:",checked,"profiles;",certified,"certified signs PASS",
      flush=True)

for d in range(3,81):
    for A in range(-3*d,3*d+1):
        lower=max(1,(A*A+d)//(d+1))
        for e in (-1,1):
            if endpoint_cert(d,A,e):
                assert lower+e*A>=0
print("integer rounding PASS",flush=True)

for d in (8,11,15,17):
    for m in (3,5):
        lo=None
        for s in (-1,1):
            bg=[(1,1,d),(1,-1,d),(2,1,d),(m,s,1)]
            A=endpoint(bg,d)
            expected=(
                (-1)**(d//4)*comb(d,d//4) if d%4==0 else
                s*(-1)**((d-m)//4)*comb(d,(d-m)//4)
                if d>=m and (d-m)%4==0 else 0
            )
            assert A==expected and (abs(A)<=1 or abs(A)>=d)
            p=pvals(bg,d)[d]
            for e in (-1,1):
                value=F(bg+[(d,e,1)],d)
                assert value==p+e*A>=0
                lo=value if lo is None else min(lo,value)
        print("endpoint family",d,m,"minimum F",lo,flush=True)

for r in (8,12,16):
    bg=[(1,1,r),(1,-1,r),(2,1,r),(3,1,1)]
    old=F(bg+[(r,1,1),(r,-1,1)],r+3)
    new=F(bg+[(3,1,1),(3,-1,1)],6)
    P=pvals(bg,r+3)
    assert old-new==sum(P[7:r+4])>=r-3 and new>0
    print("neutral family",r,"F",old,"reduced F",new,flush=True)

@lru_cache(None)
def inv(ns):
    if not ns:
        return 1
    if sum(ns)%2 or 2*max(ns)>sum(ns):
        return 0
    row={0:1}
    for n in ns:
        rr={}
        for j,a in row.items():
            for ell in range(abs(j-n),j+n+1,2):
                rr[ell]=rr.get(ell,0)+a
        row=rr
    return row.get(0,0)

def even(word):
    ans=0
    for mask in range(1<<len(word)):
        a=[]
        b=[]
        sign=1
        for j,(n,e) in enumerate(word):
            if mask>>j&1:
                a.append(n)
                sign*=e
            else:
                b.append(n)
        ans+=sign*inv(tuple(sorted(a)))*inv(tuple(sorted(b)))
    return ans

rng=Random(8811)
for j in range(20):
    bg=[(n,rng.choice((-1,1))) for n in (1,3,5)]
    l=j%3
    p=9-2*l
    n=max(3,l)+(j%(p-max(3,l)+1))
    rest=bg+[(n,1),(n,-1)]
    ep=1
    for _,e in rest:
        ep*=e
    assert even(rest+[(p,ep)])==2*F(
        [(m,e,1) for m,e in rest],n+l)

for d in range(3,10):
    bg=[(d,1),(d-1,-1),(3,1)]
    for e in (-1,1):
        rest=bg+[(d,e)]
        p=sum(n for n,_ in bg)-d
        ep=1
        for _,s in rest:
            ep*=s
        assert even(rest+[(p,ep)])==2*F(
            [(n,s,1) for n,s in rest],d)
print("34 direct invariant-multiplicity bridges PASS",flush=True)

H=difference(gram([(1,1,12),(3,-1,1)],6))
H=H.extract([0,2,4,6],[0,2,4,6])
v=S.Matrix([-46,11,-4,2])
assert (v.T*H*v)[0]==-S.Rational(23876,5)
assert F([(1,1,12),(3,-1,1)],6)==114972

bg=[(1,1,2),(3,-1,1)]
assert F(bg+[(3,1,2)],4)==13
assert F(bg+[(1,1,2)],2)==20
print("parity Gram and same-sign witnesses PASS",flush=True)

cs=[(1,1,58),(1,-1,2),(2,1,1),(11,-1,1)]
result=psd(difference(gram(cs,22)))
assert result[0] is False and result[1]==22 and result[2]<0
P=pvals(cs,30)
assert P[30]-P[29]==3966287354958387096773670983308
assert sum(n*c for n,e,c in cs)-60==13
print("distance-30 negative Gram pivot, positive scalar F PASS",
      flush=True)
print("ALL CHECKS PASS",flush=True)