# (E) residual with the chord bound min(LD, metric) on both chords of W_ij = T(j,i-1) - T(j+1,i).
# Metric (FM-MECH25): T(x,y)^2 <= K D_x D_y / ((2sqrtV - R+)(2sqrtV - R-)), h = 2x+C-N, R+ = C+h, R- = |C-h|, K = (N+2-C)^2 - h^2, V=(a+1)(e+1).
exec(open('split3.py').read())
import math
from collections import Counter
def chord(x,y,D,N,V):
    C=y-x; ld=(C+1)*math.sqrt(max(D(x)*D(y),0))
    h=2*x+C-N; Rp=C+h; Rm=abs(C-h); K=(N+2-C)**2-h*h; sV=2*math.sqrt(V)
    if Rp<sV and K>=0:
        m=math.sqrt(K*D(x)*D(y)/((sV-Rp)*(sV-Rm)))
        return min(ld,m)
    return ld
st=Counter(); byq=Counter()
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                st['pairs']+=1; q=i-j-1; s=i-j
                if abs(a-e)<=1 or min(a,e)<=2 or q==0 or (q==1 and abs(a-e) in (2,3)) or (j,i)==(N-1,N+1): st['old proved']+=1; continue
                ld=s*(math.sqrt(D(j)*D(i-1))+math.sqrt(D(j+1)*D(i)))
                if D(j)-D(i)>=ld*(1+1e-12): st['LD energy drop']+=1; continue
                b=chord(j,i-1,D,N,V)+chord(j+1,i,D,N,V)
                if D(j)-D(i)>=b*(1+1e-12): st['metric chord (new)']+=1; continue
                st['RESIDUAL']+=1; byq[min(q,6)]+=1
tot=st['pairs']
for k,v in st.most_common(): print('%-20s %8d  %.2f%%'%(k,v,100*v/tot))
print('residual by q:',sorted(byq.items()))
