# Independent check of Theorem LL4 with the kernel evaluator: random four-label words with mu1>N and mu2>N.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random
random.seed(11); n=0; bad=0
while n<250:
    r=random.randint(2,6); a=random.randint(0,10); N=a+2*r-4
    kap=[random.randint(0,3*N+8) for _ in range(4)]
    if (a+sum(kap))%2: continue
    L=sorted((k+1 for k in kap),reverse=True)
    mu1=L[3]+max(0,L[0]-L[1]-L[2]); mu2=L[0]-L[1]+L[2]-L[3]
    if not (mu1>N and mu2>N): continue
    v=phi_kernel(core(tuple(kap),()),r,a); n+=1
    if v<0: bad+=1; print('LL4 FAILS',kap,r,a,v)
print('LL4 independent checks:',n,'failures',bad)
