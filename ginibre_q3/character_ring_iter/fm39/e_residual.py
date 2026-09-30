# Residual of (E) after its proved uniform regions (RF omitted, conservative):
# strips |a-e|<=1; min(a,e)<=2; q=0; q=1 with |a-e| in {2,3}; outer endpoint (N-1,N+1); LD energy-drop region (FM-SEC25).
exec(open('split3.py').read())
import math
from collections import Counter
st=Counter(); ex=[]; byq=Counter(); bygap=Counter()
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                st['pairs']+=1; q=i-j-1; s=i-j
                if abs(a-e)<=1 or min(a,e)<=2: st['strip/min<=2']+=1; continue
                if q==0: st['q=0']+=1; continue
                if q==1 and abs(a-e) in (2,3): st['q=1 gaps 2,3']+=1; continue
                if (j,i)==(N-1,N+1): st['outer endpoint']+=1; continue
                if D(j)-D(i) >= s*(math.sqrt(D(j)*D(i-1))+math.sqrt(D(j+1)*D(i)))*(1+1e-12): st['LD energy drop']+=1; continue
                st['RESIDUAL']+=1; byq[min(q,6)]+=1; bygap[min(abs(a-e),12)]+=1
                if len(ex)<5: ex.append((a,e,j,i))
tot=st['pairs']
for k,v in st.most_common(): print('%-18s %8d  %.2f%%'%(k,v,100*v/tot))
print('residual by q (6 = >=6):',sorted(byq.items())); print('residual by |a-e| (12 = >=12):',sorted(bygap.items())); print('examples:',ex)
