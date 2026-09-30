# Which t is needed on sweep >= pi pairs (q >= 2): counts of (M_i: t = X_i^2), (M_j: t = X_j^2), (t = 2V), exact integers.
import sys
from collections import Counter
exec(open('e_sweep_split.py').read().split("AM=int(sys.argv[1])")[0])
AM=int(sys.argv[1]); st=Counter(); ex=[]
for a in range(3,AM+1):
    for e in range(3,a-1):
        c=row(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
        Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
        for j in range((N+1)//2,N+1):
            cj,Bj=C_(j),Bv[j]; imax=j
            for k in range(j+1,N+2):
                if cj*Bv[k]-Bj*C_(k)>0: imax=k
                else: break
            for i in range(imax+1,N+2):
                if i-j-1<=1: continue
                E=D[j]-D[i]; x=2*j-N; y=2*i-N
                Mt=lambda t: 0<t<4*V and 4*t*(4*V-t)*E*E>=(4*t*D[j]+(x*x-t)*Av[j]**2)*(4*t*D[i]+(y*y-t)*Av[i]**2)
                key=(Mt(y*y),Mt(x*x),Mt(2*V)); st[key]+=1
                if key==(False,False,True) and len(ex)<10: ex.append((a,e,j,i,x,y,4*V))
for k,v in sorted(st.items(),key=lambda z:-z[1]): print('(t=X_i^2, t=X_j^2, t=2V) =',k,v)
print('only-2V examples (a,e,j,i,X_j,X_i,4V):',ex)
