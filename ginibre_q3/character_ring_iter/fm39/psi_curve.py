# (E) as a polygon inequality for the curve psi_k = (c_k, B_k), B_k = c_(k-1) + c_(k+1):
#   W_ij = psi_j ^ psi_i,  D_k - D_(k+1) = psi_k ^ psi_(k+1),
#   left turn at psi_(k+1)  <=>  D_k - D_(k+2) - W_(k+2,k) >= 0  ((E) at q = 1, minus sign).
# Coverage of (E) (both signs) by: psi-sweep < pi (both signs from q=1 + OL via the convex polygon argument),
# sweep in [pi, 2pi] (minus sign trivial).  Also test an LD-type bound |W_ij| <= (i-j) sqrt((D_j-D_(j+1))(D_(i-1)-D_i)).
exec(open('split3.py').read())
import math
from collections import Counter
st=Counter(); ex=[]
def wedge(p,q): return p[0]*q[1]-p[1]*q[0]
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); B=lambda k: C_(k-1)+C_(k+1)
        psi=lambda k:(C_(k),B(k))
        for k in range((N+1)//2,N):
            assert wedge(psi(k),psi(k+1))==D(k)-D(k+1)
            assert wedge(psi(k+1)[0]-psi(k)[0:1][0] if False else (psi(k+1)[0]-psi(k)[0], psi(k+1)[1]-psi(k)[1]),(psi(k+2)[0]-psi(k+1)[0],psi(k+2)[1]-psi(k+1)[1]))==D(k)-D(k+2)-(B(k+2)*C_(k)-C_(k+2)*B(k))
        for j in range((N+1)//2,N+1):
            for i in range(j+2,N+2):
                Wij=B(i)*C_(j)-C_(i)*B(j); assert Wij==wedge(psi(j),psi(i))
                st['pairs q>=1']+=1
                vs=[psi(k) for k in range(j,i+1)]
                if any(v==(0,0) for v in vs): st['zero vector']+=1; continue
                sw=sum(math.atan2(wedge(vs[k],vs[k+1]),vs[k][0]*vs[k+1][0]+vs[k][1]*vs[k+1][1]) for k in range(len(vs)-1))
                if sw<math.pi: st['sweep<pi (both signs from q=1)']+=1
                elif sw<=2*math.pi: st['pi<=sweep<=2pi (minus sign trivial)']+=1
                else: st['sweep>2pi']+=1
                lhs=abs(Wij); g=(D(j)-D(j+1))*(D(i-1)-D(i))
                if g>=0 and lhs*lhs<=(i-j)**2*g: st['LD-psi holds']+=1
                else:
                    st['LD-psi fails']+=1
                    if len(ex)<4: ex.append((a,e,j,i,lhs,g))
tot=st['pairs q>=1']
for k,v in st.most_common(): print('%-40s %8d  %.1f%%'%(k,v,100*v/tot))
print('LD-psi failures (a,e,j,i,|W|,product):',ex)
