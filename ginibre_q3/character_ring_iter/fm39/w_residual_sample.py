# Sampled exhaustiveness check of the W criteria beyond the grid: random rows (r, a), all consumer windows of the row.
# Progress: time-based heartbeat with rows done/total, windows, residual count.
import math, random, time, sys
exec(open('w_residual4.py').read().split("from collections import Counter")[0])
random.seed(int(sys.argv[1])); rows=int(sys.argv[2])
from collections import Counter
st=Counter(); res=[]; t0=time.time(); last=t0
for rr in range(rows):
    r=random.randint(3,40); a=random.randint(0,400); e=2*r-3; N=a+e; d=a-e
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
            st['RESIDUAL']+=1
            if len(res)<10: res.append((r,a,x,Cw,round(u,2)))
    st['rows']+=1
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat rows %d/%d windows %d residual %d'%(rr+1,rows,st['windows'],st['RESIDUAL']),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed',round(time.time()-t0,1)); print('residual examples (r,a,x,C,u):',res)
