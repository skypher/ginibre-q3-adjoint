# (E) residual with every criterion: strips, min(a,e)<=2, q<=1 (EQ1), outer endpoint, LD/metric chords (both signs),
# q=2 minus sign (OL on (1-z^2)P), and the definite binary form (ratio-free, per sign).
exec(open('e_metric.py').read().split("st=Counter()")[0])
exec(open('e_binary_form.py').read().split("qmax=")[0])
from collections import Counter
st=Counter(); byq=Counter(); ex=[]
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; d=a-e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        Dd=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                st['pairs']+=1; q=i-j-1; s=i-j
                if abs(a-e)<=1 or min(a,e)<=2 or q<=1 or (j,i)==(N-1,N+1): st['strips/min/q<=1/endpoint']+=1; continue
                ld=s*(math.sqrt(Dd(j)*Dd(i-1))+math.sqrt(Dd(j+1)*Dd(i)))
                b=min(ld, chord(j,i-1,Dd,N,V)+chord(j+1,i,Dd,N,V))
                if Dd(j)-Dd(i)>=b*(1+1e-12): st['LD/metric']+=1; continue
                ok=True
                for sg in (1,-1):
                    if q==2 and sg==1: continue          # minus-W sign at q=2: OL on (1-z^2)P
                    A,B,Cq=form(d,N,j,q,sg)
                    if not(B*B-4*A*Cq<0 and A>0): ok=False
                if ok: st['definite binary form']+=1; continue
                st['RESIDUAL']+=1; byq[min(q,8)]+=1
                if len(ex)<6: ex.append((a,e,j,i))
tot=st['pairs']
for k,v in st.most_common(): print('%-28s %8d  %.3f%%'%(k,v,100*v/tot))
print('residual by q:',sorted(byq.items()),' examples:',ex)
