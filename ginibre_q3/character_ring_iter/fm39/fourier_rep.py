# Fourier representation for (anti)reciprocal rows, N even.  e even: F(t) = sum_k c_k cos(m_k t); e odd: F(t) = sum_k c_k sin(m_k t),
# m_k = k - N/2.  Claim: u = t - p, v = t + p (same formulas for both parities),
#   D_k      = sigma * mean_(t,p) F(t)F(p) cos(m_k u)(1 - cos v)                         [sigma: sign from the sine series]
#   T(x,x+C) = sigma * mean F(t)F(p) cos(M u) 2 sin((C+1)v/2) sin(v/2),  M = x + C/2 - N/2
import numpy as np
from math import comb
def row(a,e):
    c=[0]*(a+e+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return c
res={}
for (a,e) in [(4,2),(6,4),(3,5),(8,2),(5,5),(7,3),(2,8),(9,1)]:
    c=row(a,e); N=a+e
    if N%2: continue
    K=512; t=2*np.pi*np.arange(K)/K; sig=1  # same formula for both parities (derivation: c_j c_k = mean F(t)F(p) cos(m_j t - m_k p))
    base=np.cos if e%2==0 else np.sin
    F=sum(c[k]*base((k-N/2)*t) for k in range(N+1))
    # normalization: c_k recovered as scale*mean(F*base(m t))
    T1,P1=np.meshgrid(t,t,indexing='ij'); FF=np.outer(F,F); C_=lambda k: c[k] if 0<=k<=N else 0
    err=0
    for k in range(N+1):
        m=k-N/2; Dex=C_(k)**2-C_(k-1)*C_(k+1)
        Dn=sig*np.mean(FF*np.cos(m*(T1-P1))*(1-np.cos(T1+P1)))
        err=max(err,abs(Dn-Dex)/(1+abs(Dex)))
    for x in range(N+1):
        for Cw in range(1,N-x+1):
            M=x+Cw/2-N/2; Tex=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
            Tn=sig*np.mean(FF*np.cos(M*(T1-P1))*2*np.sin((Cw+1)*(T1+P1)/2)*np.sin((T1+P1)/2))
            err=max(err,abs(Tn-Tex)/(1+abs(Tex)))
    res[(a,e)]=err
print('max relative error per row:',{k:float('%.2e'%v) for k,v in res.items()})
