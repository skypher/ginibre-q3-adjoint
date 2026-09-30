# q = 2 plus sign, ratio-free: for every (d, N, j): Q_plus positive definite, or the W-form PSD (then W >= 0 at the actual row).
exec(open('e_binary_form.py').read().split("qmax=")[0])
import sys
from collections import Counter
def Wform(d,N,j,q):
    def Wv(p0,p1):
        c=window_vals(d,N,j,q,p0,p1); i=j+q+1; C=lambda k: c[k]
        return (C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
    A=Wv(1,0); Cq=Wv(0,1); B=Wv(1,1)-A-Cq
    return A,B,Cq
AM=int(sys.argv[1]); st=Counter(); ex=[]
for a in range(0,AM+1):
    for e in range(0,AM+1):
        N=a+e; d=a-e
        for j in range((N+1)//2, N-1):
            st['(d,N,j)']+=1
            A,B,Cq=form(d,N,j,2,-1)
            if B*B-4*A*Cq<0 and A>0: st['Q_plus PD']+=1; continue
            a1,b1,c1=Wform(d,N,j,2)
            if a1>=0 and c1>=0 and b1*b1-4*a1*c1<=0: st['W-form PSD']+=1; continue
            st['neither']+=1
            if len(ex)<6: ex.append((a,e,j,(A,B,Cq),(a1,b1,c1)))
print(dict(st)); print('neither examples:',[x[:3] for x in ex])
