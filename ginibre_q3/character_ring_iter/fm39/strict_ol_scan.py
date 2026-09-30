# Strict OL scan: for 2k > N, k <= N, is delta_k = D_k - D_(k+1) > 0 whenever psi_k = (c_k, B_k) != 0?  Exact integers.
import sys
from math import comb
M=int(sys.argv[1]) if len(sys.argv)>1 else 100
flat=[]; n=0; zero_psi=0
for a in range(0,M+1):
    for e in range(0,M+1):
        N=a+e; c=[0]*(N+1)
        for u in range(e+1):
            s=(-1)**u*comb(e,u)
            for k in range(a+1): c[k+u]+=s*comb(a,k)
        C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        for k in range(N//2+1,N+1):
            n+=1
            if C_(k)==0 and C_(k-1)+C_(k+1)==0: zero_psi+=1; continue
            dk=D(k)-D(k+1)
            if dk<=0: flat.append((a,e,k,dk))
print('steps with 2k>N checked:',n,'zero psi:',zero_psi,'nonpositive delta with psi != 0:',len(flat),flat[:10])
