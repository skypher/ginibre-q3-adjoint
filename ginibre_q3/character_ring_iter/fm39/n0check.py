# Independent check of Astra FM-MECH1 suggestion 1 (N0): for a core w (at most 2r general h_k, k>=2, plus hat S_p),
# P(m) = D(m) * phi_r(h_1^(2m+eps) w) / phi_r(h_1^(2m)) is a polynomial in m with nonnegative Newton coefficients.
# phi_r(w) = (1/2) E[w (x-y)^(2r)] over two independent semicircles; computed from Catalan moments (no closed form used).
import sys, itertools, time
from fractions import Fraction as Fr
from math import comb
from functools import lru_cache
def cat(n): return comb(2*n,n)//(n+1)
# polynomials in x,y as dict {(i,j):c}
def mul(A,B):
    C={}
    for (a,b),u in A.items():
        for (c,d),v in B.items():
            k=(a+c,b+d); C[k]=C.get(k,0)+u*v
    return {k:v for k,v in C.items() if v}
def add(*As):
    C={}
    for sgn,A in As:
        for k,v in A.items(): C[k]=C.get(k,0)+sgn*v
    return {k:v for k,v in C.items() if v}
ONE={(0,0):1}; X={(1,0):1}; Y={(0,1):1}
@lru_cache(None)
def U1(n,var):
    if n<0: return ()
    if n==0: return tuple(ONE.items())
    V=X if var==0 else Y
    return tuple(add((1,mul(V,dict(U1(n-1,var)))),(-1,dict(U1(n-2,var)))).items())
def U(n,var): return dict(U1(n,var))
@lru_cache(None)
def h1(k):  # h_k = sum_{p+q=k} U_p(x) U_q(y)
    A={}
    for p in range(k+1): A=add((1,A),(1,mul(U(p,0),U(k-p,1))))
    return tuple(A.items())
def h(k): return dict(h1(k))
def shat(p): return add((1,U(p,0)),(1,U(p,1)))
def E(A):  # E over semicircle x semicircle
    t=Fr(0)
    for (i,j),c in A.items():
        if i%2==0 and j%2==0: t+=c*cat(i//2)*cat(j//2)
    return t
@lru_cache(None)
def dpow(r): 
    A=ONE
    for _ in range(2*r): A=mul(A,{(1,0):1,(0,1):-1})
    return tuple(A.items())
@lru_cache(None)
def spow(n):
    A=ONE
    for _ in range(n): A=mul(A,{(1,0):1,(0,1):1})
    return tuple(A.items())
def phi(r,A,n):  # phi_r(h_1^n A)
    return E(mul(mul(A,dict(spow(n))),dict(dpow(r))))/2
def newton(vals):
    out=[]; row=list(vals)
    while row:
        out.append(row[0]); row=[row[j+1]-row[j] for j in range(len(row)-1)]
    return out
def rf(a,n):
    p=1
    for t in range(n): p*=a+t
    return p
def test(r,w,deg_w):
    # parity eps = deg_w mod 2 ; L bound: b+c <= (deg in s + eps)/2 + deg in d/2 <= (deg_w+1)//2 ... use L=deg_w//2+1 generous
    eps=deg_w%2; L=(deg_w+eps)//2
    K=2*L+3
    vals=[]
    for m in range(K+2):
        Rm=phi(r,w,2*m+eps)/phi(r,ONE,2*m)
        Dm=rf(m+r+2,L)*rf(m+r+3,L)
        vals.append(Dm*Rm)
    nw=newton(vals)
    # polynomial check: Newton coefficients beyond degree 2L must vanish
    poly_ok=all(v==0 for v in nw[2*L+1:])
    coeffs=nw[:2*L+1]
    return poly_ok, min(coeffs), coeffs, vals[0]
if __name__=='__main__':
    t0=time.time()
    # reproduce the agent's four cases
    print('agent cases:')
    for name,r,w,dg in [('h2^4',2,mul(mul(h(2),h(2)),mul(h(2),h(2))),8),('h2^6',3,mul(mul(mul(h(2),h(2)),mul(h(2),h(2))),mul(h(2),h(2))),12),
                        ('S3^2 h2^2',1,mul(mul(shat(3),shat(3)),mul(h(2),h(2))),10),('S3 S2 h3',2,mul(mul(shat(3),shat(2)),h(3)),8)]:
        ok,mn,co,v0=test(r,w,dg); print(f'  {name} r={r}: polynomial={ok} min Newton coeff={mn} (P(0)={v0})',flush=True)
