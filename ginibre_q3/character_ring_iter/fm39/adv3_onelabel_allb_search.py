"""Exact all-b profile search; no files."""
import argparse
from math import comb
argparse.ArgumentParser(description=__doc__).parse_args()

def row(a,e):
    N=a+e; c=[1]
    if N: c.append(a-e)
    for j in range(1,N):
        h=(a-e)*c[j]-(N-j+1)*c[j-1]
        assert h%(j+1)==0
        c.append(h//(j+1))
    return c

def ballot(m,n):
    if m<n or (m-n)%2: return 0
    h=(m-n)//2
    return comb(m,h)-(comb(m,h-1) if h else 0)

profiles=coefficients=residual=zeros=0
least=None
jobs=[(e,a,24) for e in range(1,13) for a in range(13)]
jobs += [(3,97,32),(4,128,32),(64,64,32),
         (31,1,32),(31,17,32),(1,127,32)]

for e,a,B in jobs:
    L=e+a; c=row(a,e)
    now=[0]*(L+1); prev=[]
    for j in range(0,L+1,2):
        now[L-j]=c[j]*comb(j,j//2)//(j//2+1)
    for b in range(B+1):
        D=L+2*b
        ns=range(D%2,D+1,2)
        us=[sum(now[m]*ballot(m,n)
                for m in range(n,D+1,2)) for n in ns]
        assert min(us)>=0,(e,a,b,min(us))
        profiles+=1; coefficients+=len(us)
        for n,val in zip(ns,us):
            if b>=3 and n>=7 and (D-n)//2>=3:
                residual+=1; zeros+=val==0
                if val>0 and (least is None or val<least[0]):
                    least=(val,e,a,b,n,(D-n)//2)
        if b==B: break
        get=lambda p,j:p[j] if 0<=j<len(p) else 0
        nxt=[]
        for j in range(D+3):
            den=D+4-j
            num=(2*(L-2-j)*get(now,j)
                 +4*b*(get(prev,j-2)+2*get(prev,j)))
            assert num%den==0
            nxt.append(get(now,j-2)+num//den)
        prev,now=now,nxt

print(profiles,coefficients,residual,zeros,least)
# 4098 86838 58654 0 (32,3,4,3,7,3)
