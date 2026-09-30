# Does P_C(x)>=0 hold for general real-rooted rows?  Test several families.
import random
from fractions import Fraction as Fr
def row(roots):
    c=[Fr(1)]
    for rho in roots:
        c=[ (c[k] if k<len(c) else 0) + rho*(c[k-1] if k>=1 else 0) for k in range(len(c)+1)]
    return c
def Pmin(c):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    T=lambda j,k: C_(j)*C_(k)-C_(j-1)*C_(k+1)
    worst=None
    for x in range(N+1):
        for Cw in range(1,N-x+1):
            P=sum(T(k,k) for k in range(x,x+Cw+1))-T(x,x+Cw)
            scale=sum(T(k,k) for k in range(x,x+Cw+1))
            if worst is None or P<worst[0]: worst=(P,x,Cw)
    return worst
random.seed(1)
fams={
 'roots +-1 only (control)':lambda: [1]*random.randint(0,10)+[-1]*random.randint(0,10),
 'roots in {+1,-lam}':lambda: (lambda l:[1]*random.randint(1,8)+[-l]*random.randint(1,8))(Fr(random.randint(1,9),random.randint(1,9))),
 'positive roots (all c>=0)':lambda: [Fr(random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(2,12))],
 'mixed-sign random roots':lambda: [Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(2,12))],
 'roots +-t pairs':lambda: sum(([t,-t] for t in [Fr(random.randint(1,5),random.randint(1,5)) for _ in range(random.randint(1,5))]),[]),
}
for name,gen in fams.items():
    neg=0; tot=0; ex=None
    for _ in range(300):
        c=row(gen()); w=Pmin(c); tot+=1
        if w and w[0]<0: neg+=1; ex=ex or (w, [str(v) for v in c[:6]])
    print(name,': rows',tot,'rows with a negative window',neg, 'example' if ex else '', ex if ex else '')
print('--- controls ---')
def quadfac(b,cc):  # multiply by z^2 + b z + cc -> as row factor (1 + b z + cc z^2)
    return (b,cc)
def rowq(roots,quads):
    c=row(roots)
    for (b,q) in quads:
        c=[ (c[k] if k<len(c) else 0) + b*(c[k-1] if 1<=k<=len(c) else 0) + q*(c[k-2] if 2<=k<=len(c)+1 else 0) for k in range(len(c)+2)]
    return c
for name,gen in {
 'one complex pair (disc<0)': lambda: ([Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(1,8))],[(Fr(random.randint(-3,3)),Fr(random.randint(3,9)))]),
 'random integer rows': None}.items():
    neg=0; tot=0; ex=None
    for _ in range(300):
        if gen is None: c=[Fr(random.randint(-9,9)) for _ in range(random.randint(3,12))]
        else:
            rt,qs=gen(); qs=[(b,q) for (b,q) in qs if b*b<4*q]; c=rowq(rt,qs)
        w=Pmin(c); tot+=1
        if w and w[0]<0: neg+=1; ex=ex or w
    print(name,': rows',tot,'with a negative window',neg,ex)
