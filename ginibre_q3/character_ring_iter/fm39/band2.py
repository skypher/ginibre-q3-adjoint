# Cancellation ratio phi / (1/2 sum |F_pq W_pq|) for balanced three-factor words, by regime of labels vs support N = a + 2r.
import itertools, math, collections
from fractions import Fraction as Fr
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
def parts(F,r,a):
    e=2*r; N=a+e; c=cvec(a,e); C=lambda k: c[k] if 0<=k<=N else 0
    tot=0; ab=0
    for (p,q),v in F.items():
        if (N+p+q)%2: continue
        P,Q=max(p,q),min(p,q); j=(N+P-Q)//2; i=(N+P+Q)//2+1
        W=(C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
        tot+=v*W; ab+=abs(v*W)
    return tot//2, ab/2
res=collections.defaultdict(list)
for r in range(2,9):
    for a in range(0,25):
        N=a+2*r
        for kap in itertools.combinations_with_replacement(range(1,19),3):
            t=sum(k%2 for k in kap)
            if (a+t)%2 or max(kap) > sum(kap)+a-max(kap): continue
            val,ab=parts(core(kap,()),r,a)
            if ab==0: continue
            lab=max(kap); reg='small' if lab*lab<=N else ('large' if lab>=N else 'band')
            res[reg].append((val/ab,kap,r,a,val))
for reg in ('small','band','large'):
    L=sorted(res[reg]); 
    if not L: continue
    qs=[L[int(q*(len(L)-1))][0] for q in (0,0.01,0.1,0.5)]
    print(f'{reg:6s}: n={len(L):6d}  cancellation ratio min={qs[0]:.2e} p1={qs[1]:.2e} p10={qs[2]:.2e} median={qs[3]:.2e}; worst {L[0][1:]}')
