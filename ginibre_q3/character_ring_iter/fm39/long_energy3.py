# (i) margin of LE on base rows; (ii) LE on general real-rooted rows; (iii) unimodality of |g_k|^2 on base rows.
exec(open('window_passes.py').read().split("stat={")[0])
exec(open('window_realrooted.py').read().split("random.seed(1)")[0])
import math, random
def LE_scan(c):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0; res=[]
    for x in range(0,N+1):
        for Cw in range(3,N-x+1):
            T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
            if T<=0: continue
            g=[(float(C_(k)),float(C_(k-1))) for k in range(x,x+Cw+2)]
            sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
            if sw<=2*math.pi: continue
            Ds=float(sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+Cw+1)))
            res.append(Ds/(math.hypot(*g[0])*math.hypot(*g[-1])))
    return res
base=[]; uni_fail=0; rows=0
for r in range(2,9):
    for a in range(0,40):
        c=list(cv(a,2*r-3)); base+=LE_scan(c); rows+=1
        N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
        rad=[C_(k)**2+C_(k-1)**2 for k in range(0,N+2)]
        # unimodal: no interior strict local minimum
        uni_fail+= any(rad[k]<rad[k-1] and rad[k]<rad[k+1] for k in range(1,N+1))
print('base rows: LE windows',len(base),' min ratio %.4f'%min(base),' rows with non-unimodal |g_k|^2:',uni_fail,'of',rows)
random.seed(4); gen=[]; nrows=0
for _ in range(400):
    c=row([Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(4,14))])
    gen+=LE_scan(c); nrows+=1
print('general real-rooted rows:',nrows,' LE windows',len(gen),' failures',sum(1 for t in gen if t<1),' min ratio', min(gen) if gen else None)
