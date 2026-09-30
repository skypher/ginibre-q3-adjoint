# Does Theorem OL (D_k >= D_(k+1) for 2k > N') hold for the rows R_m^+/- P = (1 -/+ z^m)(1+z)^a(1-z)^e, N' = N+m?
# Only the indices consumed by (E) matter: k = j + (q+m)/2 with j >= N/2, i.e. 2k >= N + q + m >= N + 2m... check all 2k > N'.
exec(open('split3.py').read())
from collections import Counter
st=Counter(); ex={}
def Dr(c,k):
    N=len(c)-1; C_=lambda t: c[t] if 0<=t<=N else 0
    return C_(k)**2-C_(k-1)*C_(k+1)
for a in range(0,31):
    for e in range(0,31):
        c=list(cv(a,e)); N=a+e
        for m in range(1,12):
            for sg in (1,-1):
                R=[0]*(m+1); R[0]=1; R[m]=-sg      # sg=+1: 1 - z^m ; sg=-1: 1 + z^m
                d=[0]*(N+m+1)
                for i,x in enumerate(c):
                    d[i]+=x; d[i+m]+=-sg*x
                Np=N+m
                for k in range(Np//2+1, Np+1):
                    if 2*k<=Np: continue
                    consumed = 2*k >= N+2*m     # indices actually used by (E): k = j+(q+m)/2 >= N/2 + m
                    key=('1-z^m' if sg==1 else '1+z^m', 'consumed' if consumed else 'other')
                    st[key]+=1
                    if Dr(d,k)<Dr(d,k+1): st[key+('FAIL',)]+=1; ex.setdefault(key,(a,e,m,k))
for k,v in sorted(st.items()): print(k,v)
print('first failures:',ex)
