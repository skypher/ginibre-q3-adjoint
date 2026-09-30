# (E) residual, sign by sign, after: strips, min(a,e)<=2, q=0, q=1 (OL on (1-/+z)P, both signs), q=2 minus-W sign (OL on (1-z^2)P),
# outer endpoint, LD energy drop and metric chords (both signs).
exec(open('e_metric.py').read().split("st=Counter()")[0])
from collections import Counter
st=Counter(); byq=Counter(); ex=[]
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                st['pairs']+=1; q=i-j-1; s=i-j
                if abs(a-e)<=1 or min(a,e)<=2 or q<=1 or (j,i)==(N-1,N+1): st['proved (incl. new q=1)']+=1; continue
                ld=s*(math.sqrt(D(j)*D(i-1))+math.sqrt(D(j+1)*D(i)))
                b=min(ld, chord(j,i-1,D,N,V)+chord(j+1,i,D,N,V))
                both = D(j)-D(i)>=b*(1+1e-12)
                if both: st['LD/metric']+=1; continue
                if q==2:
                    B=lambda k: C_(k-1)+C_(k+1); W=B(i)*C_(j)-C_(i)*B(j)
                    # minus-W sign proved by OL on (1-z^2)P; need the plus-W sign D_j - D_i + W >= 0 by the chord bound on |W|?  use: plus sign holds if W >= 0 (then D_j-D_i+W >= D_j-D_i >= 0)
                    if W>=0: st['q=2 (minus by OL, plus since W>=0)']+=1; continue
                st['RESIDUAL']+=1; byq[min(q,6)]+=1
                if len(ex)<5: ex.append((a,e,j,i))
tot=st['pairs']
for k,v in st.most_common(): print('%-40s %8d  %.2f%%'%(k,v,100*v/tot))
print('residual by q:',sorted(byq.items()),' examples:',ex)
