# For each (E) pair in the open region, check whether every term of the multiplier identity is >= 0 (termwise OL-cyc),
# for both signs.  Terms: delta_(j+q)[(z^i0 -/+ z^(q-i0)) c], i0 < q/2, plus middle term (q even, plus sign).
exec(open('split3.py').read())
from collections import Counter
def dlr(d,k):
    N=len(d)-1; C_=lambda t: d[t] if 0<=t<=N else 0
    Dk=lambda t: C_(t)**2-C_(t-1)*C_(t+1)
    return Dk(k)-Dk(k+1)
st=Counter(); ex=[]
for a in range(3,41):
    for e in range(3,41):
        if abs(a-e)<2: continue
        c=list(cv(a,e)); N=a+e
        def mult(i0,m,sg):
            d=[0]*(N+i0+m+1)
            for t,x in enumerate(c): d[t+i0]+=x; d[t+i0+m]+=-sg*x
            return d
        for j in range((N+1)//2,N+1):
            for i in range(j+2,N+2):
                q=i-j-1
                for sg in (1,-1):
                    terms=[]
                    for i0 in range(0,(q+1)//2):
                        if 2*i0>=q: break
                        terms.append(dlr(mult(i0,q-2*i0,sg),j+q))
                    if q%2==0 and sg==-1:
                        d=[0]*(N+q//2+1)
                        for t,x in enumerate(c): d[t+q//2]+=2*x
                        terms.append(dlr(d,j+q)/2)
                    st['pairs x signs']+=1
                    if all(t>=0 for t in terms): st['termwise OL-cyc']+=1
                    else:
                        st['some term negative']+=1
                        if len(ex)<5: ex.append((a,e,j,i,sg,terms))
print(dict(st)); print('examples (a,e,j,i,sign,terms):',ex)
