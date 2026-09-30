# FM3 falsification via the kernel reduction: phi_r(w h_1^a) = (1/2) sum_{p,q} F_pq(w) W(p,q; a, e=2r),
# W(p,q) = B_i c_j - c_i B_j with j=(N+p-q)/2, i=(N+p+q)/2+1 for p>=q (symmetric in p,q at even e), c = coeffs of (1+z)^a(1-z)^e.
import sys, time, json, os, itertools
from math import comb
from collections import defaultdict
def mul1(A,B):
    C=defaultdict(int)
    for (a,b),x in A.items():
        for (c,d),y in B.items():
            xy=x*y
            for e_ in range(abs(a-c),a+c+1,2):
                for f in range(abs(b-d),b+d+1,2): C[(e_,f)]+=xy
    return {k:v for k,v in C.items() if v}
def hk(k): return {(p,k-p):1 for p in range(k+1)}
def sk(p): 
    d={(p,0):1}; d[(0,p)]=d.get((0,p),0)+1; return d
def cvec(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        s=(-1)**i*comb(e,i)
        for j in range(a+1): c[i+j]+=s*comb(a,j)
    return c
def phi_kernel(F,r,a):
    e=2*r; N=a+e; c=cvec(a,e)
    C=lambda k: c[k] if 0<=k<=N else 0
    tot=0
    for (p,q),v in F.items():
        if (N+p+q)%2: continue
        P,Q=max(p,q),min(p,q); j=(N+P-Q)//2; i=(N+P+Q)//2+1
        W=(C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
        tot+=v*W
    assert tot%2==0
    return tot//2
def core(hs,ss):
    F={(0,0):1}
    for k in hs: F=mul1(F,hk(k))
    for p in ss: F=mul1(F,sk(p))
    return F
if __name__=='__main__' and sys.argv[1]=='check':
    exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
    bad=0; n=0
    for hs,ss in [((2,),()),((3,2),()),((2,2,2),()),((),(3,)),((),(2,2)),((2,),(3,)),((4,2),(2,)),((3,),(2,3))]:
        F=core(hs,ss)
        for r in range(1,5):
            for a in range(0,7):
                w=ONE
                for k in hs: w=mul(w,h(k))
                for p in ss: w=mul(w,shat(p))
                d=phi(r,w,a); kv=phi_kernel(F,r,a); n+=1
                if d!=kv: bad+=1; print('MISMATCH',hs,ss,r,a,d,kv) if bad<5 else None
    print('kernel formula vs direct phi:',n,'cases, mismatches',bad)
