# Three-factor gamma <= N branch: phi as a binary quadratic form in (c_x0, c_(x0+1)) via the recurrence; census of definiteness.
from fractions import Fraction as Fr
from collections import Counter
exec(open('split3.py').read())
def row_from(d,N,x0,p0,p1,lo,hi):
    c={x0:Fr(p0),x0+1:Fr(p1)}
    for k in range(x0+1,hi): c[k+1]=(d*c[k]-(N-k+1)*c[k-1])/(k+1)
    for k in range(x0,lo,-1): c[k-1]=(d*c[k]-(k+1)*c[k+1])/(N-k+1)
    return c
def phi_form(a,e,A,B,C):
    N=a+e; d=a-e; al=(N+A-B-C)//2; ga=al+B; tau=ga+1
    lo=min(al-A-2,-2); hi=max(tau+C+3,N+3)
    x0=max(0,min(N-1,(N)//2))
    def val(p0,p1):
        c=row_from(d,N,x0,p0,p1,lo,hi)
        C_=lambda k: c.get(k,0) if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1)
        P=lambda x: sum(D(k) for k in range(x,x+C+1))-C_(x)*C_(x+C)+C_(x-1)*C_(x+C+1)
        Ac=lambda x: C_(x)-C_(x+C); Bc=lambda x: C_(x-1)-C_(x+C+1)
        return P(al)-P(tau)-(Ac(al)*Bc(tau)-Bc(al)*Ac(tau))
    f10=val(1,0); f01=val(0,1); f11=val(1,1)
    return f10, f11-f10-f01, f01
st=Counter(); chk=0
for r in range(2,6):
    e=2*r-3
    for a in range(0,26):
        N=a+e
        for w in range(2,N+1):
            for v in range(w,N+3):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2
                    if not(al>=0 and al+B<=N): continue
                    if abs(a-e)<=1 or a==0: continue
                    Aq,Bq,Cq=phi_form(a,e,A,B,C); st['gamma<=N words (outside strips, a>0)']+=1
                    if chk<200:   # sanity: form evaluated at the true (c_x0, c_x0+1) equals phi
                        c=cv(a,e); x0=max(0,min(N-1,N//2)); p0,p1=c[x0],(c[x0+1] if x0+1<=N else 0)
                        assert Aq*p0*p0+Bq*p0*p1+Cq*p1*p1==phi3(r,a,u,v,w); chk+=1
                    if Bq*Bq-4*Aq*Cq<0 and Aq>0: st['definite']+=1
                    elif Bq*Bq-4*Aq*Cq<=0 and Aq>=0 and Cq>=0: st['semidefinite']+=1
print(dict(st),' sanity checks passed:',chk)
