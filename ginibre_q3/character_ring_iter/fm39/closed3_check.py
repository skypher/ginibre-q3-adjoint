# Independent check of the three-factor closed form phi = M - R against the kernel evaluator (fm3kern), random cases.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random
def closed(r,a,kap):
    e=2*r-3; N=a+e; c=cvec(a,e); C=lambda k: c[k] if 0<=k<=N else 0
    D=lambda k: C(k)**2-C(k-1)*C(k+1)
    A,B,Cc=sorted((k+1 for k in kap),reverse=True)
    if (N+A+B+Cc)%2: return 0
    al=(N+A-B-Cc)//2; be=al+Cc; ga=al+B; de=al+B+Cc
    Q=lambda x,y,z,t:(x+t)*(y+z)-x*t-y*z
    M=sum(D(k) for k in range(al,be+1))-sum(D(k) for k in range(ga+1,de+2))
    R=Q(C(al),C(be),C(ga),C(de+2))-Q(C(al-1),C(be+1),C(ga+1),C(de+1))
    return M-R
random.seed(7); n=0; bad=0
for _ in range(600):
    r=random.randint(2,8); a=random.randint(0,20); kap=tuple(random.randint(1,16) for _ in range(3))
    if (a+sum(kap))%2: continue
    v=phi_kernel(core(kap,()),r,a); w=closed(r,a,kap); n+=1
    if v!=w: bad+=1; print('MISMATCH',r,a,kap,v,w) if bad<4 else None
print('closed form vs kernel evaluator:',n,'random cases, mismatches',bad)
