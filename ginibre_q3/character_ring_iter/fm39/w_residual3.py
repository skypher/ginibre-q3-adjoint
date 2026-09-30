# Residual of W after all proved criteria: C<=2 certificates; T<=0; WS (sweep<=2pi); G0B window region
# 4C(1-rho_x)(1-rho_y)>=1; LD+energy (sqrt(Dmax/Dmin) not in (1,C)).  Consumer windows: e odd, 2x >= N-C, x+C <= N.
exec(open('window_passes.py').read().split("stat={")[0])
import math
from collections import Counter
st=Counter(); ex=[]
for r in range(2,11):
    e=2*r-3
    for a in range(0,60):
        N=a+e; d=a-e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+1)]
        rho=lambda k: abs(d)/(2*math.sqrt((k+1)*(N-k+1)))
        for Cw in range(3,N+1):
            for x in range(0,N-Cw+1):
                if 2*x<N-Cw: continue
                y=x+Cw; st['windows']+=1
                T=C_(x)*C_(y)-C_(x-1)*C_(y+1)
                if T<=0: st['T<=0']+=1; continue
                mx,mn=max(D[x],D[y]),min(D[x],D[y]); u=math.sqrt(mx/mn)
                if e==1: st['(ii) at e=1']+=1; continue
                if min(a,e)<=2 or abs(a-e)<=1: st['G0E+(E) proved rows']+=1; continue
                if 2*x<=N and (N-2*x+1)*math.sqrt(D[x]/D[y])>=2*x+Cw-N and D[x]>=D[y]: st['(ii) plateau']+=1; continue
                if u<=1+1e-15 or u>=Cw-1e-12: st['LD+energy']+=1; continue
                rx,ry=rho(x),rho(y)
                if rx<1 and ry<1 and 4*Cw*(1-rx)*(1-ry)>=1: st['G0B region']+=1; continue
                g=[(float(C_(k)),float(C_(k-1))) for k in range(x,y+2)]
                sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
                if sw<=2*math.pi: st['WS']+=1; continue
                # G0E telescoping: every step (y, y+Cw+1) must be a proved (E) instance (LD energy-drop region or q=1 gap 2,3 or outer endpoint)
                def eproved(j,i):
                    if i>=N+2: return True
                    if 2*j<N: j=N-j
                    if j>=i: return True
                    s_=i-j; Dj,Di,Dj1,Di1=D[j],(D[i] if i<=N else 0),(D[j+1] if j+1<=N else 0),D[i-1]
                    return Dj-Di >= s_*(math.sqrt(Dj*Di1)+math.sqrt(Dj1*Di))*(1+1e-12) or (i-j-1==1 and abs(a-e) in (2,3)) or (j,i)==(N-1,N+1) or i==j+1
                if all(eproved(yy,yy+Cw+1) for yy in range(x,N+1)): st['G0E steps in proved (E)']+=1; continue
                st['RESIDUAL']+=1
                S=sum(D[x:y+1])
                if len(ex)<6: ex.append((r,a,x,Cw,round(u,3),round(S/(Cw+1)/math.sqrt(D[x]*D[y]),4)))
tot=st['windows']
for k,v in st.most_common(): print('%-12s %8d  %.2f%%'%(k,v,100*v/tot))
print('residual examples (r,a,x,C,u,(ii)-ratio):',ex)
