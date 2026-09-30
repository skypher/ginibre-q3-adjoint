# Coverage of G0 windows by the proven criteria: WS (sweep <= 2pi or T <= 0), WL (pass condition), each for c and the twist.
import math
exec(open('window_passes.py').read().split("stat={")[0])
def criteria(cw):  # returns set of proven labels for this row window (list of float points g)
    g=cw; sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
    T=g[0][0]*g[-1][1]-g[0][1]*g[-1][0]
    if T<=0: return 'T<=0'
    if sw<math.pi: return 'convex'
    if sw>2*math.pi:
        beta=sw%(2*math.pi); thx=math.atan2(g[0][1],g[0][0]); the=math.atan2(g[-1][1],g[-1][0])
        r1=pass_radius(g,beta,the); gr=[(p[0],-p[1]) for p in g[::-1]]; rL=pass_radius(gr,beta,-thx)
        if (r1 is not None and r1>=math.hypot(*g[-1])*(1-1e-12)) or (rL is not None and rL>=math.hypot(*g[0])*(1-1e-12)): return 'pass'
    return None
from collections import Counter
cnt=Counter(); unc=[]
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        seen={}
        for w in range(1,N+2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,Cw=u+1,v+1,w+1; al=(N+A-B-Cw)//2; be=al+Cw; ga=al+B
                    if not(0<=al<be<=N and ga>=N+1): continue
                    key=(al,Cw)
                    if key not in seen:
                        g=[(float(C_(k)),float(C_(k-1))) for k in range(al,al+Cw+2)]
                        gt=[(float((-1)**k*C_(k)),float((-1)**(k-1)*C_(k-1))) for k in range(al,al+Cw+2)]
                        lab=criteria(g)
                        if lab is None:
                            if Cw%2==0:
                                lt=criteria(gt); lab=('twist-'+lt) if lt else None
                            else:
                                Tt=gt[0][0]*gt[-1][1]-gt[0][1]*gt[-1][0]
                                lab='twist-odd-T>=0' if Tt>=0 else None
                        seen[key]=lab
                    lab=seen[key]; cnt[lab or 'UNCOVERED']+=1
                    if lab is None and len(unc)<8: unc.append((r,a,(u,v,w),al,Cw))
tot=sum(cnt.values())
for k,v in cnt.most_common(): print('%-18s %8d  %.2f%%'%(k,v,100*v/tot))
print('uncovered examples (r,a,word,alpha,C):',unc)
