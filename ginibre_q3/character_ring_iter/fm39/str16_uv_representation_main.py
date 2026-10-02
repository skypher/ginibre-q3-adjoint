# FM-STR16 (main agent): the (uv, uv^{-1}) representation of FM3.
# With u Haar in SU(2) and v = diag(e^{i psi}, e^{-i psi}) such that v^2 is Haar, X = tr(uv) and Y = tr(uv^{-1}) are
# independent semicircles.  Hence Phi = sum_k ||P (x)_i e_{k_i}||^2 G_eps(k), sum k = 0, with P the projection onto
# invariants and G_eps(k) = sum_F eps^F c(sum_{i in F} k_i), c(0)=1, c(+-2)=-1/2, else 0.
# Part 1: exact identity check on all lists with labels <= 3, length 2..5.  Part 2: G takes negative values.
# Test G_eps(k) = sum_F eps^F c(sum_{i in F} k_i) >= 0, c(0)=1, c(+-2)=-1/2, else 0, for weight vectors k with sum k = 0,
# k_i in {-n_i, -n_i+2, ..., n_i}.  Also verify Phi = sum_k ||P w_k||^2 G(k) numerically via an independent formula is NOT
# done here; first just look for negative G.
import itertools
from fractions import Fraction as Fr
def c(s): return Fr(1) if s==0 else (Fr(-1,2) if abs(s)==2 else Fr(0))
def G(k,eps):
    N=len(k); tot=Fr(0)
    for F in range(1<<N):
        s=0; e=1
        for i in range(N):
            if F>>i&1: s+=k[i]; e*=eps[i]
        tot+=e*c(s)
    return tot
# Verify Phi = sum_k E_u[prod_i D^{n_i}_{k_i k_i}(u)] * G_eps(k) exactly on small lists.
import itertools
from fractions import Fraction as Fr
from math import comb, factorial
from collections import defaultdict
# polynomials in a, ab(=abar), b, bb as dict {(p,q,r,s):coef}
def pmul(P,Q):
    R=defaultdict(Fr)
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            R[tuple(x+y for x,y in zip(k1,k2))]+=v1*v2
    return dict(R)
def ppow(P,e):
    R={(0,0,0,0):Fr(1)}
    for _ in range(e): R=pmul(R,P)
    return R
def diag_coeff(n,p):
    # coefficient of z1^p z2^(n-p) in (a z1 - bb z2)^p (b z1 + ab z2)^(n-p)
    # (a z1 - bb z2)^p = sum_i C(p,i) a^i (-bb)^(p-i) z1^i z2^(p-i)
    # (b z1 + ab z2)^(n-p) = sum_j C(n-p,j) b^j ab^(n-p-j) z1^j z2^(n-p-j)
    R=defaultdict(Fr)
    for i in range(p+1):
        j=p-i
        if j<0 or j>n-p: continue
        R[(i,n-p-j,j,p-i)]+=comb(p,i)*comb(n-p,j)*(-1)**(p-i)
    return dict(R)
def haar(P):
    t=Fr(0)
    for (pa,pab,pb,pbb),v in P.items():
        if pa==pab and pb==pbb: t+=v*Fr(factorial(pa)*factorial(pb),factorial(pa+pb+1))
    return t
def phi(word):
    # direct Walsh form via tables
    d={(0,0):1}
    for z in word:
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for cc in range(abs(a-n),a+n+1,2): q[cc,b]+=v
            for cc in range(abs(b-n),b+n+1,2): q[a,cc]+=e*v
        d=q
    return d.get((0,0),0)
ok=0; bad=[]
for N in range(2,6):
    for ns in itertools.combinations_with_replacement(range(1,4),N):
        if sum(ns)%2: continue
        for eps in itertools.product((1,-1),repeat=N):
            if eps.count(-1)%2: continue
            tot=Fr(0)
            for ps in itertools.product(*[range(n+1) for n in ns]):
                k=[2*p-n for p,n in zip(ps,ns)]
                if sum(k)!=0: continue
                P={(0,0,0,0):Fr(1)}
                for n,p in zip(ns,ps): P=pmul(P,diag_coeff(n,p))
                w=haar(P)
                if w: tot+=w*G(k,eps)
            word=[e*n for e,n in zip(eps,ns)]
            if tot==phi(word): ok+=1
            else: bad.append((word,tot,phi(word)))
print("identity checks ok",ok,"mismatches",len(bad),bad[:3])

worst=(Fr(10),None); cnt=0; neg=0
for N in range(2,7):
    for ns in itertools.combinations_with_replacement(range(1,5),N):
        if sum(ns)%2: continue
        for eps in itertools.product((1,-1),repeat=N):
            if eps.count(-1)%2: continue
            for k in itertools.product(*[range(-n,n+1,2) for n in ns]):
                if sum(k)!=0: continue
                gg=G(k,eps); cnt+=1
                if gg<0: neg+=1
                if gg<worst[0]: worst=(gg,(ns,eps,k))
print("weight vectors",cnt,"negative G",neg,"min",worst)
assert not bad and ok==240
print("FM-STR16 IDENTITY PASS; G is sign-indefinite")
