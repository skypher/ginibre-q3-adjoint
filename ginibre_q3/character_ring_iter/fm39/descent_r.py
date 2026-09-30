import sys
exec(open('em2r.py').read().split("R=int(sys.argv[1])")[0])
X={(1,0):1,(0,1):-1}
def parts(n,m,l):
    if n==0: yield (); return
    if l==0: return
    for k in range(min(n,m),0,-1):
        for r in parts(n-k,k,l-1): yield (k,)+r
cacheP={}
def prod(kap):
    kap=tuple(sorted([k for k in kap if k],reverse=True))
    if kap in cacheP: return cacheP[kap]
    P={(0,0):1} if not kap else mul1(prod(kap[:-1]),hk(kap[-1]))
    cacheP[kap]=P; return P
def rules(lam):
    s=len(lam); out={}
    L=list(lam)
    def e(i,j):
        k=L[:]; k[i]-=1; k[j]-=1; return tuple(k)
    if s>=2:
        out['edge smallest two']=e(s-2,s-1)
        out['edge largest two']=e(0,1)
        out['edge largest-smallest']=e(0,s-1)
        out['merge smallest two']=tuple(L[:-2])+(L[-2]+L[-1],)
    if L[0]>=2: out['largest -2']=tuple([L[0]-2]+L[1:])
    if L[-1]>=2: out['smallest -2']=tuple(L[:-1]+[L[-1]-2])
    if s>=2: out['drop smallest two if equal'] = tuple(L[:-2]) if L[-1]==L[-2] else None
    return {k:v for k,v in out.items() if v is not None}
N=int(sys.argv[1])
for R in (2,3,4,5):
    Q={(0,0):1}
    for _ in range(2*R): Q=mul1(Q,X)
    phi=lambda k: sum(v*Q.get(key,0) for key,v in prod(k).items())//2
    stats={}
    for n in range(2,N+1,2):
        for lam in parts(n,n,99):
            v=phi(lam)
            for name,k2 in rules(lam).items():
                d=v-phi(k2); st=stats.setdefault(name,[0,0,None])
                st[0]+=1
                if d<0:
                    st[1]+=1
                    if st[2] is None: st[2]=(lam,d)
    print(f"r={R}:",{k:(f"{v[1]}/{v[0]} fail",v[2]) for k,v in stats.items()})
