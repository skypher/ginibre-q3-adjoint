# Four-factor support region (outer pairing 12|34): s* > N and L1+L2-L3-L4 > N  =>  corr = 0 and
# phi = sum_(d,d') phi(hatS_d hatS_d'), each a two-label T+W value.  Independent check with phi_row.
exec(open('general_row_words.py').read().split("random.seed(")[0])
from collections import Counter
def cg_(p,q): return range(abs(p-q),p+q+1,2)
def lam(I):
    s=sum(I); lo=max(0,2*max(I)-s)
    return lo if (lo-s)%2==0 else lo+1
st=Counter()
import itertools
for r in range(2,6):
    e=2*r-4
    for a in range(0,9):
        c=mkrow(a,e,[]); N=len(c)-1
        for Ls in itertools.combinations_with_replacement(range(3,12),4):
            L=sorted(Ls,reverse=True)
            if (N+sum(L))%2: continue
            sstar=min(L[i]+lam([L[j] for j in range(4) if j!=i]) for i in range(4))
            if not(sstar>N and L[0]+L[1]-L[2]-L[3]>N): continue
            phi=phi_row(c,L,[],e%2)
            main=sum(phi_row(c,[],[d,dp],e%2) for d in cg_(L[0],L[1]) for dp in cg_(L[2],L[3]))
            mu1=L[3]+max(0,L[0]-L[1]-L[2]); mu2=L[0]-L[1]+L[2]-L[3]
            st['region words']+=1; st['corr=0']+= main==phi; st['phi>=0']+= phi>=0
            st['outside LL4']+= not(mu1>N and mu2>N)
print(dict(st))
