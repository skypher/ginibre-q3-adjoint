# Test the three-factor formula, OL and (E) on general (anti)reciprocal real-rooted rows
# P(z) = (1+z)^a (1-z)^e prod_i (1 + t_i z + z^2), |t_i| >= 2  (roots rho, 1/rho real).
import random
from fractions import Fraction as Fr
from math import comb
def mulrow(c,f):
    out=[Fr(0)]*(len(c)+len(f)-1)
    for i,x in enumerate(c):
        for j,y in enumerate(f): out[i+j]+=x*y
    return out
def mkrow(a,e,ts):
    c=[Fr(1)]
    for _ in range(a): c=mulrow(c,[1,1])
    for _ in range(e): c=mulrow(c,[1,-1])
    for t in ts: c=mulrow(c,[1,t,1])
    return c
def tools(c):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
    return N,C_,D
def Phi3(c,al,B,Cw):   # three-factor closed form (1): P_C(l)-P_C(tau)-(A_l B_tau - B_l A_tau)
    N,C_,D=tools(c)
    P=lambda x: sum(D(k) for k in range(x,x+Cw+1))-C_(x)*C_(x+Cw)+C_(x-1)*C_(x+Cw+1)
    Ac=lambda x: C_(x)-C_(x+Cw); Bc=lambda x: C_(x-1)-C_(x+Cw+1)
    tau=al+B+1
    return P(al)-P(tau)-(Ac(al)*Bc(tau)-Bc(al)*Ac(tau))
random.seed(11)
def randt(): 
    s=random.choice([-1,1]); return s*(2+Fr(random.randint(0,12),random.randint(1,6)))
stats={}
def rec(k,bad,ex=None):
    s=stats.setdefault(k,[0,0,None]); s[0]+=1
    if bad: s[1]+=1; s[2]=s[2] or ex
# sanity: with ts=[] and e odd, Phi3 must equal the direct word value
exec(open('split3.py').read())
for r in range(2,5):
    for a in range(0,6):
        e=2*r-3; c=mkrow(a,e,[]); N=a+e
        for w in range(1,5):
            for v in range(w,7):
                for u in range(v,9):
                    if (a+u+v+w)%2: continue
                    A,B,Cw=u+1,v+1,w+1; al=(N+A-B-Cw)//2
                    assert Phi3(c,al,B,Cw)==phi3(r,a,u,v,w),(r,a,u,v,w)
print('closed form (1) == direct phi on base rows: ok')
for trial in range(1500):
    e=random.choice([1,3,5]); a=random.randint(0,6); ts=[randt() for _ in range(random.randint(1,3))]
    c=mkrow(a,e,ts); N=len(c)-1
    # three-factor: all label triples realizable at this degree, including words with gamma<=N
    for Cw in range(2,N+2):
        for B in range(Cw,N+3):
            for A in range(B,B+Cw+1):
                if (N+A-B-Cw)%2: continue
                al=(N+A-B-Cw)//2
                v=Phi3(c,al,B,Cw)
                rec('three-factor, anti-reciprocal real-rooted', v<0, (a,e,[str(t) for t in ts],A,B,Cw,str(v)))
    # OL-type on this row: D_k >= D_{k+1} for 2k>N
    N_,C_,D=tools(c)
    for k in range(N//2+1,N+1): rec('OL (D_k>=D_k+1, 2k>N)', D(k)<D(k+1), (a,e,[str(t) for t in ts],k))
for trial in range(800):   # two-label (E) on reciprocal rows (e even)
    e=random.choice([0,2,4]); a=random.randint(0,6); ts=[randt() for _ in range(random.randint(1,3))]
    c=mkrow(a,e,ts); N,C_,D=tools(c)
    for p in range(0,N+1):
        for q in range(0,p+1):
            if (N+p+q)%2 or p+q>N: continue
            j=(N+p-q)//2; i=(N+p+q)//2+1
            W=(C_(i-1)+C_(i+1))*C_(j)-C_(i)*(C_(j-1)+C_(j+1))
            rec('(E) D_j-D_i>=|W|, reciprocal real-rooted', D(j)-D(i)<abs(W), (a,e,[str(t) for t in ts],p,q))
    for k in range(N//2+1,N+1): rec('OL on reciprocal rows', D(k)<D(k+1), (a,e,[str(t) for t in ts],k))
for k,v in stats.items(): print(k,': tested',v[0],'failures',v[1],'first',v[2])
