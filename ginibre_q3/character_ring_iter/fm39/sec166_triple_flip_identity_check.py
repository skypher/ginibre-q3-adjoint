# Check: for B, p, and a same-parity pair (a,b) in B, C = B - a - b, Lambda = B + (sigma p):
#   g_p(B) = T1 + (D_ab + D_ap + D_bp)/2,  T1 = sum_{c in CG(a,b)} sum_{s in CG(p,c)} g_s(C),
# where D_uv = eps_v A_{n_u n_v}(Lambda - u - v) are the flip quarter-differences.
import random
from collections import defaultdict
def cg(a,b): return range(abs(a-b),a+b+1,2)
def table(L):
    A={(0,0):1}
    for z in L:
        n=abs(z); e=1 if z>0 else -1; T=defaultdict(int)
        for (s,t),v in A.items():
            for c in cg(s,n): T[c,t]+=v
            for c in cg(t,n): T[s,c]+=e*v
        A={k:v for k,v in T.items() if v}
    return A
def D(Lam,u,v):
    C=list(Lam); C.remove(u); C.remove(v); a=table(C).get((abs(u),abs(v)),0); return a if v>0 else -a
rng=random.Random(7); checked=0
for _ in range(400):
    B=[rng.choice([-1,1])*rng.randint(1,7) for _ in range(rng.randint(3,8))]
    # pair-free: one sign per label
    sg={}; B=[sg.setdefault(abs(z),1 if z>0 else -1)*abs(z) for z in B]
    W=sum(map(abs,B)); p=rng.randint(max(map(abs,B)),W)
    if (W-p)%2: continue
    sigma=(-1)**sum(z<0 for z in B); Lam=B+[sigma*p]
    idx=[(i,j) for i in range(len(B)) for j in range(i+1,len(B)) if (B[i]+B[j])%2==0]
    if not idx: continue
    i,j=rng.choice(idx); a,b=B[i],B[j]; C=[z for k,z in enumerate(B) if k not in (i,j)]
    TC=table(C); gB=table(B).get((p,0),0)
    T1=sum(TC.get((s,0),0) for c in cg(abs(a),abs(b)) for s in cg(p,c))
    rhs2=2*T1+D(Lam,a,b)+D(Lam,a,sigma*p)+D(Lam,b,sigma*p)
    assert 2*gB==rhs2,(B,p,a,b,gB,T1)
    checked+=1
print("identity g_p(B) = T1 + (D_ab + D_ap + D_bp)/2 checked on",checked,"random cases")
