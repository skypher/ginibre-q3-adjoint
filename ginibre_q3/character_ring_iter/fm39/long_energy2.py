# Long sweeps (> 2 pi, T > 0), C >= 3: sum D >= |g_x||g_end| (=> sum D >= T), for c or (even C) its twist.
exec(open('window_passes.py').read().split("stat={")[0])
import math
from collections import Counter
st=Counter(); worst=None; fails=[]
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        for x in range(0,N+1):
            for Cw in range(3,N-x+1):
                T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
                if T<=0: continue
                g=[(float(C_(k)),float(C_(k-1))) for k in range(x,x+Cw+2)]
                sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
                if sw<=2*math.pi: continue
                st['long']+=1
                Ds=sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+Cw+1))
                R=math.hypot(*g[0])*math.hypot(*g[-1])
                ok = Ds>=R
                st['energy ok']+=ok
                if not ok:
                    fails.append((round(Ds/R,3),r,a,x,Cw,round(sw/(2*math.pi),2)))
fails.sort()
print(dict(st)); print('failures (ratio,r,a,x,C,turns):',len(fails),fails[:6])
print('failures by a:',Counter(f[2] for f in fails).most_common(8))
