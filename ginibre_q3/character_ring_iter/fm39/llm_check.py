# Independent check of Theorem LLm (m = 5, 6) with the kernel evaluator.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random, itertools
def lam(Ls):
    T=sum(Ls); M=max(Ls); d=max(0,2*M-T)
    return d if (d-T)%2==0 else d+1
random.seed(5); n=0; bad=0; tries=0
while n<120 and tries<200000:
    tries+=1
    m=random.choice([5,6]); r=random.randint((m+1)//2,5); a=random.randint(0,6); N=a+2*r-m
    kap=[random.randint(0,4*N+10) for _ in range(m)]
    if (a+sum(kap))%2: continue
    L=[k+1 for k in kap]; ok=True
    for size in range(1,m//2+1):
        for S in itertools.combinations(range(m),size):
            Sc=[i for i in range(m) if i not in S]
            if lam([L[i] for i in S])+lam([L[i] for i in Sc])<=N: ok=False; break
        if not ok: break
    if not ok: continue
    v=phi_kernel(core(tuple(kap),()),r,a); n+=1
    if v<0: bad+=1; print('LLm FAILS',kap,r,a,v)
print('LLm independent checks (m=5,6):',n,'failures',bad)
