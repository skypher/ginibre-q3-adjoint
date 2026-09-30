# Explore four-factor words at large labels: phi vs one-label main term tau(U_ABCD), and the sign of the 2|2 pairing terms.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import itertools, random
def cvW(r,a,e):
    N=a+e; c=cvec(a,e); C=lambda k: c[k] if 0<=k<=N else 0
    def W(p,q):
        if (N+p+q)%2 or p+q>N: return 0
        P,Q=max(p,q),min(p,q); i=(N+P+Q)//2+1; j=(N+P-Q)//2
        return (C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
    return W
def cg(p,q): return range(abs(p-q),p+q+1,2)
def prodU(labels):
    d={0:1}
    for L in labels:
        nd={}
        for l,m in d.items():
            for k in cg(l,L): nd[k]=nd.get(k,0)+m
        d=nd
    return d
def WW(W,f,g): return sum(m*n*W(p,q) for p,m in f.items() for q,n in g.items())
random.seed(3); stats={'phi<0':0,'pair_sum<0':0,'n':0}; worst=None
for _ in range(300):
    r=random.randint(2,5); a=random.randint(0,10); e=2*r-4; N=a+e
    kap=tuple(sorted(random.randint(max(1,N-1),N+6) for _ in range(4)))
    if (a+sum(kap))%2: continue
    A=[k+1 for k in kap]; W=cvW(r,a,e)
    tau=WW(W,prodU(A),{0:1})
    pairs=0
    for (i,j),(k,l) in [((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]:
        pairs+=WW(W,prodU([A[i],A[j]]),prodU([A[k],A[l]]))
    three=sum(WW(W,prodU([A[x] for x in S]),{A[y]:1}) for y in range(4) for S in [tuple(z for z in range(4) if z!=y)])
    phi=phi_kernel(core(kap,()),r,a)
    pred=tau-three+pairs
    assert phi==pred, (kap,r,a,phi,pred)
    stats['n']+=1; stats['phi<0']+=phi<0; stats['pair_sum<0']+=pairs<0
    if pairs<0 and (worst is None or pairs/tau<worst[0]): worst=(pairs/tau,kap,r,a,tau,pairs,three)
print(stats); print('most negative pairing sum relative to tau:',worst)
