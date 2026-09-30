# Independent check of Theorem LL with the kernel evaluator (fm3kern): identity (2) and positivity on large-label triples.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import itertools, random
def phi_word(kap,r,a): return phi_kernel(core(kap,()),r,a)
def mu_list(A,B,C):
    mu={}
    for d in range(abs(A-B),A+B+1,2):
        for l in range(abs(d-C),d+C+1,2): mu[l]=mu.get(l,0)+1
    return mu
random.seed(1); bad=0; n=0
for _ in range(400):
    r=random.randint(2,7); a=random.randint(0,12); N=a+2*r-3
    lo=max(1,N-1); kap=tuple(sorted(random.randint(lo,lo+8) for _ in range(3)))
    if (a+sum(kap))%2: continue
    A,B,C=[k+1 for k in kap]
    main=sum(m*phi_word((l-1,),r-1,a) for l,m in mu_list(A,B,C).items() if 1<=l<=N)
    edge=(C==N and A==B)+(B==N and A==C)+(A==N and B==C)
    v=phi_word(kap,r,a); n+=1
    if v!=main+edge or v<0: bad+=1; print('FAIL',kap,r,a,v,main,edge)
print('Theorem LL identity (2) independent checks:',n,'failures',bad)
