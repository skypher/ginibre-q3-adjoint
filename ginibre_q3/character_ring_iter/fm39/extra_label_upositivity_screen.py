# U-positivity of g_{k,a,b}(x) = E_y[(x-y)^k (x+y)^a (x^2+y^2-2)^b], y semicircle on [-2,2]
import sys
from math import comb
from functools import lru_cache
def cat(j): return comb(2*j,j)//(j+1)
def polymul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a:
            for j,b in enumerate(q): r[i+j]+=a*b
    return r
# bivariate poly as dict {(i,j):c}
def bmul(P,Q):
    R={}
    for (i,j),a in P.items():
        for (k,l),b in Q.items():
            R[(i+k,j+l)]=R.get((i+k,j+l),0)+a*b
    return R
def bpow(P,n):
    R={(0,0):1}
    for _ in range(n): R=bmul(R,P)
    return R
D={(1,0):1,(0,1):-1}; S={(1,0):1,(0,1):1}; Zp={(2,0):1,(0,2):1,(0,0):-2}
def g(k,a,b):
    F=bmul(bmul(bpow(D,k),bpow(S,a)),bpow(Zp,b))
    deg=max(i for (i,j) in F)
    p=[0]*(deg+1)
    for (i,j),c in F.items():
        if j%2==0: p[i]+=c*cat(j//2)
    return p
def ucoef(p):
    # expand monomial poly in U_n: x^m = sum_i ballot(m,i) U_{m-2i}, ballot = C(m,i)-C(m,i-1)
    out=[0]*len(p)
    for m,c in enumerate(p):
        if c:
            for i in range(m//2+1):
                out[m-2*i]+=c*(comb(m,i)-(comb(m,i-1) if i else 0))
    return out
bad=0; tot=0; monopos=0; zeros=0
K=int(sys.argv[1])
for k in range(1,K+1):
  for a in range(0,K+1):
    for b in range(0,K+1):
      if k+a+2*b>2*K+4: continue
      p=g(k,a,b); u=ucoef(p); tot+=1
      # only parity-consistent n matter: U_n coefficient with n = k+a mod 2
      if min(u)<0:
          bad+=1; print("NEG",k,a,b,u)
      if min(p)>=0: monopos+=1
print("cases",tot,"negative",bad,"monomial-positive",monopos)
print("examples"); 
for (k,a,b) in [(2,0,1),(4,0,0),(2,2,1),(3,1,2),(6,0,0),(2,0,3)]:
    print((k,a,b),"mono",g(k,a,b),"U",ucoef(g(k,a,b)))
