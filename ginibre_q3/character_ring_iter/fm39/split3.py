# Direct three-factor evaluator from the definition:
# phi_r(h_u h_v h_w h1^a) = (1/2) E[(x-y)^e (x+y)^a prod_i (U_Ai(x)-U_Ai(y))], e=2r-3,
# expanded over splits, W(p,q)=sum_k c_k m(N-k,p) m(k,q), m(i,p)=E[x^i U_p(x)] (ballot numbers).
from math import comb
from functools import lru_cache
@lru_cache(None)
def _mom(i,p):
    if i<p or (i-p)%2: return 0
    return (p+1)*comb(i,(i-p)//2)//((i+p)//2+1)
@lru_cache(None)
def cv(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return tuple(c)
@lru_cache(None)
def W(a,e,p,q):
    N=a+e
    if p+q>N: return 0
    c=cv(a,e); return sum(c[k]*_mom(N-k,p)*_mom(k,q) for k in range(q,N-p+1))
def cg(p,q): return range(abs(p-q),p+q+1,2)
def phi3(r,a,u,v,w):
    e=2*r-3; A,B,C=u+1,v+1,w+1
    tot=0
    for l1 in cg(A,B):
        for l in cg(l1,C): tot+=W(a,e,l,0)          # all in x
        tot-=W(a,e,l1,C)                            # {A,B}|{C}
    for l1 in cg(A,C): tot-=W(a,e,l1,B)
    for l1 in cg(B,C): tot-=W(a,e,l1,A)
    # remaining four splits are the swaps; (x-y)^e odd => W(Q,P)=-W(P,Q), doubling
    assert tot%1==0
    return tot   # = (1/2)*2*tot
