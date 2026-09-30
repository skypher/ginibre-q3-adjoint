# Main agent: exact half-angle determinant formula for [U_n] g_(k,a,b),
# g_(k,a,b)(x) = E_y[(x-y)^k (x+y)^a (x^2+y^2-2)^b] (one-extra-label sector).
# Checks [U_n] g = sum_i C(b,i) 2^(b-i) (u_n v_(n+2) - u_(n+2) v_n) for all k+a+2b <= K
# (argument K, default 22) and counts negative single-i terms.  Run from fm39/.
import sys
from math import comb
exec(open('extra_label_upositivity_screen.py').read().split("bad=0")[0])
def lmul(A,B):
    R={}
    for i,a in A.items():
        for j,b in B.items(): R[i+j]=R.get(i+j,0)+a*b
    return {k:v for k,v in R.items() if v}
def lpow(A,n):
    R={0:1}
    for _ in range(n): R=lmul(R,A)
    return R
tm={1:1,-1:-1}; tp={1:1,-1:1}; t2={2:1,-2:1}
def Ai(k,a,i): return lmul(lmul(lpow(tm,k),lpow(tp,a)),lpow(t2,i))
def Phi(A,n):
    u=lambda j:A.get(j,0); v=lambda j:A.get(j-2,0)+A.get(j+2,0)
    return u(n)*v(n+2)-u(n+2)*v(n)
K=int(sys.argv[1]) if len(sys.argv)>1 else 22; tot=0; bad_term=0; mism=0; worst=[]
for k in range(1,K+1):
  for a in range(0,K+1):
    for b in range(0,K//2+1):
      if k+a+2*b>K: continue
      uc=ucoef(g(k,a,b))
      for n in range(0,k+a+2*b+1):
        if (n-k-a)%2: continue
        terms=[comb(b,i)*2**(b-i)*Phi(Ai(k,a,i),n) for i in range(b+1)]
        val=uc[n] if n<len(uc) else 0
        tot+=1
        if sum(terms)!=val: mism+=1; print("MISMATCH",k,a,b,n,sum(terms),val)
        for i in range(b+1):
            ph=Phi(Ai(k,a,i),n)
            if ph<0: bad_term+=1; worst.append((ph,k,a,i,n))
print("entries",tot,"identity mismatches",mism,"negative single-i terms",bad_term)
worst.sort(); print(worst[:8])
