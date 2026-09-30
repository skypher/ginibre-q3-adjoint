# (E) <= (i)-psi and (ii)-psi, with delta_k = D_k - D_(k+1) >= 0 (Theorem OL, k >= N/2):
#   (i)-psi : |W_ij| <= (i-j) sqrt(delta_j delta_(i-1))
#   (ii)-psi: sum_(k=j..i-1) delta_k >= (i-j) sqrt(delta_j delta_(i-1))
exec(open('split3.py').read())
from collections import Counter
st=Counter(); ex=[]
for a in range(0,41):
    for e in range(0,41):
        c=cv(a,e); N=a+e; C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); B=lambda k: C_(k-1)+C_(k+1)
        dl=lambda k: D(k)-D(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+2,N+2):
                W=B(i)*C_(j)-C_(i)*B(j); g=dl(j)*dl(i-1); n=i-j; S=D(j)-D(i)
                open_region = min(a,e)>=3 and abs(a-e)>=2
                st['pairs']+=1; st['open pairs']+=open_region
                i_ok = g>=0 and W*W<=n*n*g; ii_ok = S>=0 and S*S>=n*n*g
                st['(i)-psi']+=i_ok; st['(ii)-psi']+=ii_ok; st['both']+=i_ok and ii_ok
                if open_region: st['open: both']+=i_ok and ii_ok
                if open_region and not ii_ok and len(ex)<4: ex.append((a,e,j,i,S,g))
for k,v in st.items(): print(k,v)
print('open-region (ii)-psi failures:',ex)
