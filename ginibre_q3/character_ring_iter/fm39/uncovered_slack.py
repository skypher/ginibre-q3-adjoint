# Slack sum D / T on the G0 windows left uncovered by window_cover.py's criteria, and on all long-sweep T>0 windows.
exec(open('window_cover.py').read().split("from collections import Counter")[0])
import math
rat_unc=[]; rat_long=[]
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        for x in range(0,N+1):
            for Cw in range(2,N-x+1):
                g=[(float(C_(k)),float(C_(k-1))) for k in range(x,x+Cw+2)]
                T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
                if T<=0: continue
                sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
                if sw<=2*math.pi: continue
                Ds=sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+Cw+1))
                rat=Ds/T; rat_long.append(rat)
                lab=criteria(g)
                if lab is None:
                    gt=[(float((-1)**k*C_(k)),float((-1)**(k-1)*C_(k-1))) for k in range(x,x+Cw+2)]
                    if Cw%2==0: lab=criteria(gt)
                    else:
                        Tt=gt[0][0]*gt[-1][1]-gt[0][1]*gt[-1][0]; lab='t' if Tt>=0 else None
                if lab is None: rat_unc.append((rat,r,a,x,Cw,round(sw/(2*math.pi),2)))
rat_unc.sort(); rat_long.sort()
print('long-sweep T>0 windows:',len(rat_long),' min sumD/T = %.3f'%rat_long[0],' 1st percentile %.3f'%rat_long[len(rat_long)//100])
print('uncovered windows:',len(rat_unc),' min sumD/T:',rat_unc[:5])
