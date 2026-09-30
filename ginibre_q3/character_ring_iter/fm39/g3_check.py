# Independent check of Lemma 2 identity (1) under (G0) and of positivity in (GC)/(GO), via the kernel evaluator.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import random
def cv(a,e):
    c=cvec(a,e); N=a+e; return (lambda k: c[k] if 0<=k<=N else 0)
random.seed(21); n=0; bad=0; pos_checked=0
while n<300:
    r=random.randint(2,6); a=random.randint(0,40); e=2*r-3; N=a+e
    w=random.choice([k for k in range(1,15) if (k+1)%2==0]); v=random.randint(w,40); u=random.randint(v,60)
    if (a+u+v+w)%2: continue
    A,B,C=u+1,v+1,w+1
    al=(N+A-B-C)//2; be=al+C; ga=al+B
    if not (0<=al<be<=N and ga>=N+1): continue
    c=cv(a,e); g=lambda k: c(k+1)-c(k-1); h=C//2
    formula=sum(g(al+2*i-1)*g(al+2*j-1)-g(al+2*i-2)*g(al+2*j) for i in range(1,h+1) for j in range(i,h+1))
    val=phi_kernel(core((u,v,w),()),r,a); n+=1
    if val!=formula: bad+=1; print('IDENTITY FAIL',r,a,(u,v,w),val,formula)
    m=min(a,e); d=A-B-C; Y=A-B+C
    if (m%2==1 and (m+1)*Y*Y<=2*(N-m+3)) or (d>=0 and d*d>=4*m*(N-m+3)):
        pos_checked+=1
        if val<=0: bad+=1; print('POSITIVITY FAIL',r,a,(u,v,w),val)
print('G3 identity checks:',n,' region positivity checks:',pos_checked,' failures:',bad)
