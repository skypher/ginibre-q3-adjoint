# Long sweeps (> 2 pi, T > 0): is sum D >= (|g_x|^2 + |g_end|^2)/2 (which implies sum D >= T)?
exec(open('window_passes.py').read().split("stat={")[0])
import math
worst=None; cnt=0; bad=0
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        for x in range(0,N+1):
            for Cw in range(2,N-x+1):
                T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
                if T<=0: continue
                g=[(float(C_(k)),float(C_(k-1))) for k in range(x,x+Cw+2)]
                sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
                if sw<=2*math.pi: continue
                Ds=sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+Cw+1))
                E=(C_(x)**2+C_(x-1)**2+C_(x+Cw+1)**2+C_(x+Cw)**2)/2
                cnt+=1; bad+= Ds<E
                rr=Ds/E
                if worst is None or rr<worst[0]: worst=(rr,r,a,x,Cw,round(sw/(2*math.pi),2))
print('long-sweep T>0 windows',cnt,' failing sumD >= (|g_x|^2+|g_end|^2)/2:',bad,' worst ratio',worst)
