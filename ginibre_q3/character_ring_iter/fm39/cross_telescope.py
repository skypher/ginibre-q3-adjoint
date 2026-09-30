# Cross terms telescope: W(U_X U_Z, U_Y) = sum_(d in CG(X,Z)) W(d,Y) = S(j0,i0) - S(j1+1,i1+1),
# S(p,q) = c_p c_(q-1) - c_(p-1) c_q, j(d) = (N+d-Y)/2, i(d) = (N+d+Y)/2+1, d0 = |X-Z|, d1 = X+Z.
exec(open('split3.py').read())
import random
random.seed(5); n=ok=0
for _ in range(4000):
    r=random.randint(2,6); a=random.randint(0,20); e=2*r-3; N=a+e
    X,Y,Z=[random.randint(1,14) for _ in range(3)]
    if (N+X+Y+Z)%2: continue
    c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
    S=lambda p,q: C_(p)*C_(q-1)-C_(p-1)*C_(q)
    j=lambda d:(N+d-Y)//2; i=lambda d:(N+d+Y)//2+1
    d0,d1=abs(X-Z),X+Z
    lhs=sum(W(a,e,d,Y) if d>=Y else (-1)**e*W(a,e,Y,d) for d in cg(X,Z))
    # W(p,q) for p<q via kernel symmetry W(p,q) = (-1)^e W(q,p)
    lhs2=sum(W(a,e,d,Y) for d in cg(X,Z))
    rhs=S(j(d0),i(d0))-S(j(d1)+1,i(d1)+1)
    n+=1; ok+= (lhs2==rhs)
print('cross-term telescoping checks',n,'matches',ok)
