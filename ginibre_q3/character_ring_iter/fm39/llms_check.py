# Independent check of Theorem LLm-S (words with h and hat S factors) with the kernel evaluator.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random, itertools
def lam(Ls):
    T=sum(Ls); M=max(Ls); d=max(0,2*M-T)
    return d if (d-T)%2==0 else d+1
random.seed(9); n=0; bad=0; tries=0
while n<150 and tries<300000:
    tries+=1
    m=random.randint(1,4); ns=random.randint(1,2); r=random.randint(max(1,(m+1)//2),5); a=random.randint(0,6); N=a+2*r-m
    hs=[random.randint(0,3*N+8) for _ in range(m)]; ss=[random.randint(2,3*N+8) for _ in range(ns)]
    if (a+sum(hs)+sum(ss))%2: continue
    labels=[k+1 for k in hs]+ss; J=len(labels); ok=True
    for size in range(1,J):
        for S in itertools.combinations(range(J),size):
            Sc=[i for i in range(J) if i not in S]
            if lam([labels[i] for i in S])+lam([labels[i] for i in Sc])<=N: ok=False; break
        if not ok: break
    if not ok: continue
    v=phi_kernel(core(tuple(hs),tuple(ss)),r,a); n+=1
    if v<0: bad+=1; print('LLm-S FAILS',hs,ss,r,a,v)
print('LLm-S independent checks:',n,'failures',bad)
