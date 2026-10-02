from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb
from random import Random
from collections import defaultdict
import argparse
ap=argparse.ArgumentParser(description="Fresh exact FM-MECH166 checks.")
ap.add_argument("--quick",action="store_true",help="reduce only random screens")
args=ap.parse_args()
rng=Random(107)

def cg(a,b):return range(abs(a-b),a+b+1,2)
def mu3(a,b,c,t):
    if (a+b-c-t)%2:return 0
    lo=max(abs(a-b),abs(c-t));hi=min(a+b,c+t)
    return max(0,(hi-lo)//2+1)
def lowc(m,r):return Q(1)+min(Q(r),Q(m,2))
def check_lemma1():
    checks=0;top=13 if args.quick else 18
    for a,b,c in combinations_with_replacement(range(1,top+1),3):
        for r in range(a+1):
            cval=lowc(a,r)
            for t in range(a+b+c+1):
                base=mu3(a,b,c,t)
                after=sum(mu3(a,b,c,s) for s in cg(t,2*r))
                assert Q(after)>=cval*base
                checks+=1
    for _ in range(500 if args.quick else 2000):
        a=rng.randrange(1,41);b=rng.randrange(a,1_000_001);c=rng.randrange(a,1_000_001)
        u=abs(a-b)+2*rng.randrange(a+1)
        t=abs(u-c)+2*rng.randrange(min(u,c)+1);r=rng.randrange(a+1)
        assert Q(sum(mu3(a,b,c,s) for s in cg(t,2*r)))>=lowc(a,r)*mu3(a,b,c,t)
    print("Lemma1 coefficientwise:",checks,"small and",
          500 if args.quick else 2000,"large random checks PASS")
def capacities(m):
    c=[lowc(m,r) for r in range(m+1)]
    return sum(c,Q(0)),sum((z*z for z in c),Q(0))
def alpha_closed(m):
    q=m//2
    return Q((q+1)*(3*q+(5 if m%2 else 2)),2)
def beta_closed(m):
    q=m//2
    return Q((q+1)*(2*q+3)*(8*q+13),12) if m%2 else Q((q+1)*(8*q*q+13*q+6),6)
def nu(L):
    if L==6:return 10
    return max(2**(L-2)-t-t*(L-t) for t in range(2,L+1,2))
def rho(L,m,b=None):
    if b is None:b=m
    aa,bb=capacities(m)
    if L==6:return 1-Q(10,1)/bb-(1+Q(3,m+1))/(b+1)
    J=L-2;K=2**(J-1)-1-J-comb(J,2)
    return 1-Q(nu(L),1)/bb-(1+Q(comb(J,2),1)/aa+Q(K,1)/bb)/(b+1)
def check_caps_thresholds():
    expected={6:(3,4),7:(4,4),8:(6,6),9:(8,8),10:(10,11),
      11:(13,14),12:(17,18),13:(22,22),14:(28,28),15:(35,36),16:(45,45)}
    for m in range(1,101):
        aa,bb=capacities(m)
        assert aa==alpha_closed(m) and bb==beta_closed(m)
        assert aa>=Q(3*(m+1)**2,8) and bb>=Q((m+1)**3,6)
    got={}
    for L in range(6,17):
        mf=1
        while beta_closed(mf)<nu(L):mf+=1
        mr=1
        while rho(L,mr)<=0:mr+=1
        got[L]=(mf,mr);assert got[L]==expected[L]
    print("Lemma2 closed forms and all L=6..16 thresholds PASS:",got)
def mu_table(ns):
    L=len(ns);N=1<<L;W=[0]*N;size=[0]*N;mx=[0]*N
    for s in range(1,N):
        bit=s&-s;i=bit.bit_length()-1
        W[s]=W[s^bit]+ns[i];size[s]=size[s^bit]+1;mx[s]=max(mx[s^bit],ns[i])
    m=[0]*N;m[0]=1
    for s in range(1,N):
        k=size[s];w=W[s]
        if k==1 or w%2 or 2*mx[s]>w:continue
        sub=s;v=0
        while True:
            q=w//2-W[sub]-size[sub]
            if q>=0:v+=(-1 if size[sub]&1 else 1)*comb(q+k-2,k-2)
            if not sub:break
            sub=(sub-1)&s
        assert v>=0;m[s]=v
    return m
def walsh(a):
    a=a[:];h=1
    while h<len(a):
        for i in range(0,len(a),2*h):
            for j in range(i,i+h):
                x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
        h*=2
    return a
def sign_masks(ns):
    classes=[];i=0
    while i<len(ns):
        j=i+1
        while j<len(ns) and ns[j]==ns[i]:j+=1
        classes.append(tuple(range(i,j)));i=j
    for cm in range(1<<len(classes)):
        mask=neg=0
        for q,cl in enumerate(classes):
            if cm>>q&1:
                neg+=len(cl)
                for i in cl:mask|=1<<i
        if neg%2==0:yield mask
def exhaustive(L,lo,hi):
    multisets=signings=0
    for ns in combinations_with_replacement(range(lo,hi+1),L):
        m=mu_table(ns);F=(1<<L)-1
        values=walsh([m[s]*m[F^s] for s in range(F+1)])
        for mask in sign_masks(ns):
            assert values[mask]>=0,(ns,mask,values[mask])
            signings+=1
        multisets+=1
    print(f"FM3 exhaustive L={L} labels=[{lo},{hi}] multisets={multisets} signings={signings} PASS")
def coeff2(factors,X,Y):
    rem=sum(abs(z) for z in factors);d={(0,0):1}
    for z in factors:
        n=abs(z);eps=1 if z>0 else -1;rem-=n;out=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n):
                if abs(c-X)+abs(b-Y)<=rem:out[c,b]+=v
            for c in cg(b,n):
                if abs(a-X)+abs(c-Y)<=rem:out[a,c]+=eps*v
        d={k:v for k,v in out.items() if v}
    return d.get((X,Y),0)
def signs(ns):
    cls=[];i=0
    while i<len(ns):
        j=i+1
        while j<len(ns) and ns[j]==ns[i]:j+=1
        cls.append(list(range(i,j)));i=j
    neg=[False]*len(ns)
    for cl in cls:
        if rng.randrange(2):
            for i in cl:neg[i]=True
    if sum(neg)%2:
        cl=next(c for c in cls if len(c)%2)
        for i in cl:neg[i]=not neg[i]
    return [-n if neg[i] else n for i,n in enumerate(ns)]
def random_large():
    rem={9:8,10:11,11:13,12:17,13:22,14:28,15:35,16:45}
    words=topchecks=0
    for L in range(9,17):
        ns=sorted(rem[L]+1+rng.randrange(16) for _ in range(L))
        signed=signs(ns)
        phi=coeff2(signed,0,0)
        gp=coeff2(signed[:-1],ns[-1],0)
        assert phi>=0 and phi%2==0 and phi==2*gp
        words+=1
        if L in (9,13,14,16):
            pairs=[(i,j) for i in range(L-1) for j in range(i+1,L-1)
                   if (ns[i]-ns[j])%2==0]
            i,j=max(pairs,key=lambda ij:(ns[ij[0]]+ns[ij[1]],
                                         max(ns[ij[0]],ns[ij[1]])))
            child=[z for k,z in enumerate(signed[:-1]) if k not in (i,j)]
            assert gp>=coeff2(child,ns[-1],0)
            topchecks+=1
    print("random above-threshold FM3 rows",words,
          "direct TopPair rows",topchecks,"PASS")

check_lemma1()
check_caps_thresholds()
exhaustive(9,8,12)
exhaustive(10,10,13)
random_large()
from fractions import Fraction as Q
from collections import defaultdict
from math import comb
from random import Random
import argparse
ap=argparse.ArgumentParser(description="Fresh exact kernel check of FM-MECH165 Theorem 5.")
args=ap.parse_args()
rng=Random(165)

def cg(a,b):return range(abs(a-b),a+b+1,2)

def coeff2(factors,X,Y):
    rem=sum(abs(z) for z in factors);d={(0,0):1}
    for z in factors:
        n=abs(z);eps=1 if z>0 else -1;rem-=n;out=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n):
                if abs(c-X)+abs(b-Y)<=rem:out[c,b]+=v
            for c in cg(b,n):
                if abs(a-X)+abs(c-Y)<=rem:out[a,c]+=eps*v
        d={k:v for k,v in out.items() if v}
    return d.get((X,Y),0)

def f_term(n,eps,cap):
    out=defaultdict(int)
    for h in range(min(n,cap)+1):out[h,h,0]+=1
    for q in range(n//2+1):
        d=n-2*q;coef=eps*(-1)**q*comb(n-q,q)
        for h in range(d+1):
            i=q+h;j=q+d-h
            if i<=cap and j<=cap:out[i,j,d]+=coef*comb(d,h)
    return [(i,j,k,v) for (i,j,k),v in out.items() if v]

def kernel(C,cap):
    P={(0,0,0):1}
    for z in C:
        terms=f_term(abs(z),1 if z>0 else -1,cap);N=defaultdict(int)
        for (i,j,d),v in P.items():
            for a,b,k,w in terms:
                if i+a<=cap and j+b<=cap:N[i+a,j+b,d+k]+=v*w
        P={key:v for key,v in N.items() if v}
    return P

class KT:
    def __init__(self,C,cap):
        self.cap=cap;self.poly=kernel(C,cap);self.coef=defaultdict(lambda:defaultdict(int))
        for (i,j,d),v in self.poly.items():self.coef[i,j][d]+=v
    def poly_at(self,i,j):
        if i<0 or j<0 or i>self.cap or j>self.cap:return {}
        return self.coef.get((i,j),{})
    def M(self,i,j):return sum(self.poly_at(i,j).values())
    def H(self,i,j):
        if min(i,j)<0:return 0
        if i<j:i,j=j,i
        return self.M(i,j)-self.M(i+1,j-1)
    def P(self,j):return self.H(j,j)
    def J(self,k,i,j):
        if min(i,j)<0:return Q(0)
        return sum((Q(2*v,k+d+2) for d,v in self.poly_at(i,j).items()),Q(0))

def check_lemma3():
    count=0;minF=None
    for m in (6,7):
      for x in range(m+1):
       for y in range(m-x+1):
        c1=comb(x,2)
        for eps in (-1,1):
         P1=m+c1;beta=eps*y+c1
         P2=comb(m+1,2)-x+c1*(m-2+eps*y)+comb(y,2)+2*comb(x,4)
         for tau in (-1,1):
          F=P2-2*P1-2*tau*beta
          minF=F if minF is None else min(minF,F)
          count+=1
    assert count==256 and minF==0
    increments=0
    for m in range(7,101):
     for x in range(m+1):
      c1=comb(x,2)
      for y in range(m-x+1):
       z=m-x-y
       for eps in (-1,1):
        for tau in (-1,1):
         P1=m+c1;beta=eps*y+c1
         P2=comb(m+1,2)-x+c1*(m-2+eps*y)+comb(y,2)+2*comb(x,4)
         F=P2-2*P1-2*tau*beta
         assert F>=0
         # Exact increments under appending a factor of label >=3, 2, or 1.
         def FF(mm,xx,yy,ee,tt):
          cc=comb(xx,2);bb=ee*yy+cc
          pp=comb(mm+1,2)-xx+cc*(mm-2+ee*yy)+comb(yy,2)+2*comb(xx,4)
          return pp-2*(mm+cc)-2*tt*bb
         assert FF(m+1,x,y,eps,tau)-F==m-1+c1
         inc2=FF(m+1,x,y+1,eps,tau)-F
         assert inc2==m-1+y+c1*(1+eps)-2*tau*eps
         inc1=FF(m+1,x+1,y,eps,tau)-F
         exact1=m-2+c1+x*(m-3-2*tau+eps*y)+2*comb(x,3)
         low1=m-2+c1+x*(x+z-5)+2*comb(x,3)
         assert inc1==exact1 and inc1>=low1>=0
         increments+=3
    print("Lemma3 corners",count,"m=7..100 exact recurrence increments",increments,"PASS")

def random_band(branch):
    delta=rng.randrange(8,16)
    if branch==0:a=b=delta
    elif branch==1:a=b=delta-1
    else:a,b=delta,delta-2
    cap=b
    small=[rng.choice((1,2,3,4)) for _ in range(2)]
    low=max(3,delta//2)
    high=min(cap,low+3)
    C=small+[rng.randint(low,high) for _ in range(4+rng.randrange(3))]
    B=[a,b]+C;W=sum(B);p=W-2*delta
    assert p>max(B) and max(B)<=delta and p>=6
    eps={n:(-1 if rng.randrange(2) else 1) for n in set(B)}
    SB=[eps[n]*n for n in B]
    ep=-1 if sum(z<0 for z in SB)%2 else 1
    full=SB+[ep*p]
    assert sum(z<0 for z in full)%2==0
    assert all(-z not in full for z in full)
    return delta,p,SB,full

def check_band():
    rows=Ecount=Jbounds=0
    for branch in range(3):
     for rep in range(4):
        delta,p,B,full=random_band(branch)
        L=len(full);W=sum(map(abs,B))
        assert L>=9 and delta==(W-p)//2 and p>=max(map(abs,B))
        ps=[(i,j) for i in range(len(B)) for j in range(i+1,len(B))
            if (abs(B[i])-abs(B[j]))%2==0]
        i,j=max(ps,key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                                  max(abs(B[ij[0]]),abs(B[ij[1]]))))
        a,b=sorted((abs(B[i]),abs(B[j])),reverse=True)
        expected=((delta,delta) if branch==0 else
                  (delta-1,delta-1) if branch==1 else (delta,delta-2))
        assert (a,b)==expected,(branch,delta,a,b,expected)
        ea=1 if B[i]>0 else -1;eb=1 if B[j]>0 else -1
        C=[v for k,v in enumerate(B) if k not in (i,j)]
        m=len(C);r=delta-a;s=delta-b;ell=delta-a-b-1
        assert (r,s)==((0,0) if branch==0 else (1,1) if branch==1 else (0,2))
        K=KT(C,delta+2)
        S=sum(K.P(q) for q in range(s,delta+1))-sum(K.P(q) for q in range(ell,r))
        X=eb*(K.H(delta,s)-K.H(r-1,ell))
        Y=ea*(K.H(delta,r)-K.H(s-1,ell))
        Z=ea*eb*(K.M(s,r)-K.M(delta+1,ell)-K.M(s-1,r-1)+K.M(delta,ell-1))
        child=coeff2(C,p,0)
        gp=coeff2(B,p,0)
        parent=S+X+Y+Z
        upper=eb*K.H(delta,s)+ea*K.H(delta,r)
        Q0=parent-child-upper
        assert gp==parent and ell<0 and Q0>0
        assert coeff2(full,0,0)==2*gp
        E=[]
        for k in range(2*b+1):
            v=(b+1)**2*K.J(2*b-k,s,s)+(a+1)**2*K.J(2*a-k,r,r)
            v+=2*ea*eb*(a+1)*(b+1)*K.J(a+b-k,r,s)
            e=K.J(k,delta,delta)*v
            assert e>=upper*upper
            E.append(e);Ecount+=1
        assert Q(Q0*Q0)>=E[0]

        # Lemma 3 formulas match the first two diagonal kernel rows.
        x=sum(abs(z)==1 for z in C);y=sum(abs(z)==2 for z in C)
        eps2=next((1 if z>0 else -1 for z in C if abs(z)==2),1)
        c1=comb(x,2);beta=eps2*y+c1
        P1=m+c1
        P2=comb(m+1,2)-x+c1*(m-2+eps2*y)+comb(y,2)+2*comb(x,4)
        assert K.P(1)==P1 and K.P(2)==P2
        assert abs(beta)<=P1 and P2>=2*(P1+abs(beta))

        # Lemma 4: exact polynomial in t=c^2 and its Bernstein coefficients.
        poly22=K.poly_at(2,2)
        v0=poly22.get(0,0);v1=poly22.get(2,0);v2=poly22.get(4,0)
        assert all(d in (0,2,4) for d in poly22)
        beta0=Q(v0);beta1=Q(v0)+Q(v1,2);beta2=Q(v0+v1+v2)
        assert beta0>=0 and beta1>=0 and beta2>=0
        assert beta0+beta1+beta2==3*P2
        def epslab(q):
            return next((1 if z>0 else -1 for z in C if abs(z)==q),1)
        m3=sum(abs(z)==3 for z in C);m4=sum(abs(z)==4 for z in C)
        beta1_formula=(comb(m+1,2)-x+c1*(m-2)-comb(y,2)
                       -eps2*c1*y-2*epslab(1)*epslab(3)*x*m3
                       -2*epslab(4)*m4)
        zsmall=m-x-y
        beta1_lower=(comb(m+1,2)-x+c1*(x+zsmall-2)-comb(y,2)
                     -2*max(x,1)*zsmall)
        assert beta1==beta1_formula and beta1>=beta1_lower>=0
        for bb in range(delta+3):
            assert K.J(2*bb,2,2)<=Q(3*P2,bb+3)
            Jbounds+=1

        # Theorem 5's branch-specific bounds.
        pd=K.P(delta)
        if branch==0:
            assert ea==eb
            assert Q0>=pd+1+(delta-1)*m
            assert E[0]==4*(delta+1)*pd
        elif branch==1:
            assert ea==eb
            assert Q0>=pd+(delta-2)*P2
            assert E[0]<=2*delta*pd*P2
            assert Q(Q0*Q0)>=4*(delta-2)*pd*P2>=E[0]
        else:
            assert Q0>=pd+Q(2*delta-5,2)*P2+1
            Ebound=pd*((Q(3*(delta-1)**2,delta+1)+Q(delta-1,2))*P2+delta+1)
            assert E[0]<=Ebound
            assert Q(Q0*Q0)-E[0]>=pd*Q(5*delta**2+2*delta-147,delta+1)>0
        assert gp>=child
        rows+=1
    print("Theorem5 random band rows",rows,"E_k checks",Ecount,
          "Lemma4 Bernstein bounds",Jbounds,"PASS")

check_lemma3()
check_band()
