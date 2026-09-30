# Independent check of decomposition (1) (general, all corrections) and positivity in (GX)/(GOdd) via the kernel evaluator.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random
random.seed(33); n=0; bad=0; gx=0; godd=0
while n<400:
    r=random.randint(2,6); a=random.randint(0,40); e=2*r-3; N=a+e
    w=random.randint(1,14); v=random.randint(w,30); u=random.randint(v,50)
    if (a+u+v+w)%2: continue
    A,B,C=u+1,v+1,w+1
    c=cvec(a,e); cc=lambda k: c[k] if 0<=k<=N else 0
    D=lambda k: cc(k)**2-cc(k-1)*cc(k+1)
    al=(N+A-B-C)//2; ga=al+B
    def P(x): return sum(D(k) for k in range(x,x+C+1))-cc(x)*cc(x+C)+cc(x-1)*cc(x+C+1)
    Acal=lambda x: cc(x)-cc(x+C); Bcal=lambda x: cc(x-1)-cc(x+C+1)
    tau=ga+1
    formula=P(al)-P(tau)-(Acal(al)*Bcal(tau)-Bcal(al)*Acal(tau))
    val=phi_kernel(core((u,v,w),()),r,a); n+=1
    if val!=formula: bad+=1; print('DECOMP FAIL',r,a,(u,v,w),val,formula)
    be=al+C; m=min(a,e); Y=A-B+C
    if C%2==0 and m%2==0 and m>=2 and 0<=al<be<=N and ga>=N+1 and 3*m*Y*Y<=N-m+4:
        gx+=1; bad+= val<=0
    if C%2==1 and 0<=al<be<=N and ga>=N+1 and (e+1)*(Y+1)**2<=2*(a+2):
        godd+=1; bad+= val<=0
print('decomposition (1) checks:',n,' GX region cases:',gx,' GOdd region cases:',godd,' failures:',bad)
