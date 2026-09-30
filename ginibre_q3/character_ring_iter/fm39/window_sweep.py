# Theorem WS check: for real-rooted rows, windows whose phase curve g_k=(c_k,c_{k-1}) sweeps <= 2 pi about O have P_C(x) >= 0
# (proof: Newton for P and (1-z)P + convex polygon).  Coverage of the G0 branch by sweep <= 2 pi.
import math
from math import comb
from fractions import Fraction as Fr
exec(open('split3.py').read())
def sweep_and_P(c,x,Cw):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    g=[(C_(k),C_(k-1)) for k in range(x,x+Cw+2)]
    sw=0.0; Dsum=0
    for k in range(len(g)-1):
        D=g[k][0]*g[k+1][1]-g[k][1]*g[k+1][0]; dot=g[k][0]*g[k+1][0]+g[k][1]*g[k+1][1]
        assert D>=0; Dsum+=D; sw+=math.atan2(D,dot)
    T=g[0][0]*g[-1][1]-g[0][1]*g[-1][0]
    return sw, Dsum-T
# hypotheses on base rows: D of c and of f=(1-z)c nonnegative
for a in range(0,25):
    for e in range(0,25):
        c=list(cv(a,e)); f=[(c[k] if k<len(c) else 0)-(c[k-1] if k>=1 else 0) for k in range(len(c)+1)]
        for row in (c,f):
            n=len(row)-1
            assert all(row[k]**2-(row[k-1] if k>=1 else 0)*(row[k+1] if k<n else 0)>=0 for k in range(n+1))
print('Newton hypotheses D(c)>=0, D((1-z)c)>=0 hold on all base rows a,e<25')
tot=short=0; bad_short=0; g3like=0; g3like_short=0
maxsw={}
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); m=min(a,e)
        for w in range(1,N+2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2; be=al+C; ga=al+B
                    if not(0<=al<be<=N and ga>=N+1): continue
                    sw,P=sweep_and_P(c,al,C); tot+=1
                    assert P==phi3(r,a,u,v,w) if tot%97==0 else True
                    if sw<=2*math.pi+1e-12:
                        short+=1; bad_short+= P<0
print('G0 words:',tot,' with sweep <= 2pi (covered by Theorem WS):',short,'(%.1f%%)'%(100*short/tot),' negative among covered:',bad_short)
