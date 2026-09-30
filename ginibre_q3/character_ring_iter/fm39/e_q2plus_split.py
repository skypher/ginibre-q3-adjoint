# q = 2 plus sign: is every pair covered by "definite binary form" or "W >= 0" (trivial)?  All a, e in range, j >= N/2.
exec(open('e_binary_form.py').read().split("qmax=")[0])
exec(open('split3.py').read())
import sys
from collections import Counter
AM=int(sys.argv[1]); st=Counter(); ex=[]
for a in range(0,AM+1):
    for e in range(0,AM+1):
        N=a+e; d=a-e; c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
        for j in range((N+1)//2, N-1):
            i=j+3; W=(C_(i-1)+C_(i+1))*C_(j)-C_(i)*(C_(j-1)+C_(j+1))
            st['pairs']+=1
            if W>=0: st['W>=0 (trivial)']+=1; continue
            A,B,Cq=form(d,N,j,2,-1)
            if B*B-4*A*Cq<0 and A>0: st['definite']+=1; continue
            st['neither']+=1
            if len(ex)<8: ex.append((a,e,j))
print(dict(st),'examples:',ex)
