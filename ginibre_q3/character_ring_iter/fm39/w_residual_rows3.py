# Sampled exhaustiveness check of the W criteria beyond the grid: random rows (r, a), all consumer windows of the row.
# Progress: time-based heartbeat with rows done/total, windows, residual count.
import math, random, time, sys
exec(open('w_residual4.py').read().split("from collections import Counter")[0])
exec(open('window_passes.py').read().split("stat={")[0])
def pass_ok(g,sw):
    beta=sw%(2*math.pi); thx=math.atan2(g[0][1],g[0][0]); the=math.atan2(g[-1][1],g[-1][0])
    r1=pass_radius(g,beta,the); gr=[(p[0],-p[1]) for p in g[::-1]]; rL=pass_radius(gr,beta,-thx)
    return (r1 is not None and r1>=math.hypot(*g[-1])*(1-1e-12)) or (rL is not None and rL>=math.hypot(*g[0])*(1-1e-12))
R1,R2,A1,A2=[int(v) for v in sys.argv[1:5]]; ROWS=[(r,a) for r in range(R1,R2+1) for a in range(A1,A2+1)]; rows=len(ROWS)
from collections import Counter
st=Counter(); res=[]; t0=time.time(); last=t0
for rr in range(rows):
    r,a=ROWS[rr]; e=2*r-3; N=a+e; d=a-e
    c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
    D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+1)]
    rho=lambda k: abs(d)/(2*math.sqrt((k+1)*(N-k+1)))
    V=(a+1)*(e+1); sV=2*math.sqrt(V)
    for Cw in range(3,N+1):
        for x in range(0,N-Cw+1):
            if 2*x<N-Cw: continue
            y=x+Cw; st['windows']+=1
            T=C_(x)*C_(y)-C_(x-1)*C_(y+1)
            if T<=0: continue
            mx,mn=max(D[x],D[y]),min(D[x],D[y]); u=math.sqrt(mx/mn)
            if u<=1+1e-15 or u>=Cw-1e-12: continue
            if e==1 or 2*x+Cw==N or min(a,e)<=2 or abs(a-e)<=1: continue
            h=2*x+Cw-N; n=Cw+1; Lh=max(1,n-h); eta=n*n if 2*Lh>=n else 4*Lh*(n-Lh)
            Rp=Cw+h; Rm=abs(Cw-h); Kq=(N+2-Cw)**2-h*h
            if Rp<sV and eta*(sV-Rp)*(sV-Rm)>=Kq*(1+1e-12): continue
            if 2*x<=N and (N-2*x+1)*math.sqrt(D[x]/D[y])>=2*x+Cw-N and D[x]>=D[y]: continue
            rx,ry=rho(x),rho(y)
            if rx<1 and ry<1 and 4*Cw*(1-rx)*(1-ry)>=1: continue
            g=[(float(C_(k)),float(C_(k-1))) for k in range(x,y+2)]
            sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
            if sw<=2*math.pi: continue
            if T>0 and pass_ok(g,sw): st['pass']+=1; continue
            # twisted row c_k -> (-1)^k c_k: D invariant, P_C invariant for even C; for odd C, P_C = sum D + T~
            gt=[(float((-1)**k*C_(k)),float((-1)**(k-1)*C_(k-1))) for k in range(x,y+2)]
            swt=sum(ang(gt[k],gt[k+1]) for k in range(len(gt)-1))
            if Cw%2==0 and swt<=2*math.pi: st['twist WS']+=1; continue
            Ttw=gt[0][0]*gt[-1][1]-gt[0][1]*gt[-1][0]
            if Cw%2==0 and swt>2*math.pi and Ttw>0 and pass_ok(gt,swt): st['twist pass']+=1; continue
            Tt=gt[0][0]*gt[-1][1]-gt[0][1]*gt[-1][0]
            if Cw%2==1 and Tt>=0: st['twist odd T~>=0']+=1; continue
            st['RESIDUAL']+=1
            res.append((r,a,x,Cw,round(u,2)))
    st['rows']+=1
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat rows %d/%d windows %d residual %d'%(rr+1,rows,st['windows'],st['RESIDUAL']),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed',round(time.time()-t0,1)); print('residual count',len(res)); print('by (r,a):',sorted(set((t[0],t[1]) for t in res))[:40]); print('first:',res[:8])
