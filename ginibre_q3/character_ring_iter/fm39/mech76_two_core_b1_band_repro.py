"""FM-MECH76: exact combined-row identities and uniform positivity regions."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from collections import defaultdict
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-balanced",type=int,default=160)
args=ap.parse_args()

# Original energy, forced energy, and b=1 coordinate change.
sig,X,de,p,q,f=S.symbols("sig X de p q f")
tt=(sig+X)/2
uu=(sig-X)/2
pn=(de*p-uu*q)/tt
H=lambda a,b:sig*(a*a+b*b)-2*de*a*b
ell=sig*q-de*p
assert S.cancel(H(p,q)-H(pn,p)-X*ell**2/tt**2)==0
forced=pn+f/tt
assert S.cancel(H(p,q)-H(forced,p)
 -(X*ell**2+2*uu*ell*f-sig*f*f)/tt**2)==0
D=p*p-q*pn
assert S.cancel(2*tt*D-H(p,q)-X*(p*p-q*q))==0
omega2=sig*sig-de*de
assert S.expand(H(p,q)**2-omega2*(p*p-q*q)**2
 -(de*(p*p+q*q)-2*sig*p*q)**2)==0
pm2=(de*q-(tt-1)*p)/(uu+1)
pp2=(de*pn-(uu-1)*p)/(tt+1)
U=(tt*tt+uu*uu+de*de-2)/((tt+1)*(uu+1))
V=de*(sig+1)/((tt+1)*(uu+1))
assert S.cancel(-pm2-pp2-U*p+V*(q+pn))==0

# Two-step norm certificate and exact constants.
rr,cc=S.symbols("rr cc")
aa=1-rr*rr
lam=rr*rr+cc*aa
trace=2*rr*rr+cc*cc*aa*aa
assert S.expand(lam*lam-trace*lam+rr**4
 -cc*cc*(1-cc)*aa**3)==0
assert F(13,25)**18<=F(1,65536)
assert F(13,25)**6<=F(1,36)
assert F(4,135)-F(1,60)-F(1,184320)>F(1,80)
assert F(13,25)**20<=F(1,512**2)
assert F(4,135)-F(1,60)-F(1,20480)>F(1,80)

# Universal insertion and determinant identities.
qm2,qm1,q0,qp1,qp2=S.symbols("z0:5")
lhs=2*(q0*q0-qm1*qp1)+(qm1+qp1)**2-(qm2+q0)*(q0+qp2)
rhs=qm1**2-qm2*q0+qp1**2-q0*qp2+q0**2-qm2*qp2
assert S.expand(lhs-rhs)==0
a0,a1,a2,b0,b1,b2=S.symbols("a0 a1 a2 b0 b1 b2")
assert S.expand(2*(a0*b1-a1*b0)+(a1*b2-a2*b1)
 -((2*a0-a2)*b1-a1*(2*b0-b2)))==0
print("symbolic identities and uniform constants PASS",flush=True)

def row(N,de):
    c=[1]
    old=0
    for k in range(N):
        z,rem=divmod(de*c[-1]-(N-k+1)*old,k+1)
        assert rem==0
        old=c[-1]
        c.append(z)
    return c

def data(N,de,B):
    c=row(N,de)
    def at(k):
        return c[k] if 0<=k<=N else 0
    @lru_cache(None)
    def v(h,k):
        if h<0:
            return 0
        return sum(comb(h,u)*at(k+h-2*u) for u in range(h+1))
    @lru_cache(None)
    def det(h,k):
        return v(h,k)**2-v(h,k-1)*v(h,k+1)
    @lru_cache(None)
    def turb(B0,k):
        return sum(comb(B0,h)*2**(B0-h)*det(h,k)
                   for h in range(B0+1))
    def T(k):
        return turb(B,k)
    def W(j,i):
        return sum(comb(B,h)*2**(B-h)*
                   (v(h,j)*v(h+1,i)-v(h+1,j)*v(h,i))
                   for h in range(B+1))
    def H0(k):
        return (N+2)*(at(k)**2+at(k-1)**2)-2*de*at(k)*at(k-1)
    return at,v,det,turb,T,W,H0

# Coupled recurrence and positive insertion.
recurrence_checks=insertion_checks=0
for N in range(21):
    for e in range(N+1):
        de=N-2*e
        at,v,det,turb,T,W,H0=data(N,de,4)
        for h in range(5):
            for k in range(-h,N+h+2):
                assert ((k+h+1)*v(h,k+1)
                        ==de*v(h,k)-(N+h-k+1)*v(h,k-1)
                          +4*h*v(h-1,k))
                recurrence_checks+=1
        for B in range(1,5):
            for k in range(-4,N+5):
                extra=sum(comb(B-1,h)*2**(B-1-h)*
                          (v(h,k)**2-v(h,k-2)*v(h,k+2))
                          for h in range(B))
                assert extra>=0
                assert turb(B,k)==turb(B-1,k-1)+turb(B-1,k+1)+extra
                assert turb(B,k)>=0
                insertion_checks+=1
print("coupled recurrence checks",recurrence_checks,
      "positive insertions",insertion_checks,flush=True)

# Balanced b=1 rows: all right-half intervals.
balanced=0
for L in range(1,args.max_balanced+1):
    N=2*L
    at,v,det,turb,T,W,H0=data(N,0,1)
    A=[-at(k-2)-at(k+2) for k in range(N+4)]
    B=[at(k-1)+at(k+1) for k in range(N+4)]
    ts=[T(k) for k in range(N+4)]
    def C(r):
        return comb(L,r) if 0<=r<=L else 0
    for r in range(L+1):
        assert A[2*r]==(-1)**r*(C(r-1)+C(r+1))
        assert B[2*r]==0
        assert ts[2*r]==(
            C(r)**2+C(r)*(C(r-1)+C(r+1))-C(r-1)*C(r+1))
        assert A[2*r+1]==0
        assert B[2*r+1]==(-1)**r*(C(r)-C(r+1))
        assert ts[2*r+1]==C(r)**2+C(r+1)**2
    for j in range(L,N+2):
        assert ts[j]>=ts[j+1]>=0
        amp=abs(A[j])+abs(B[j])
        amp2=abs(A[j+2])+abs(B[j+2])
        if amp and amp2:
            assert ts[j]*amp2>=ts[j+2]*amp
        for i in range(j+1,N+2):
            cross=A[j]*B[i]-B[j]*A[i]
            assert cross==W(j,i)
            assert ts[j]-ts[i]>=abs(cross)
            balanced+=1
print("balanced intervals",balanced,flush=True)

# Actual residual obstructions.
at,v,det,turb,T,W,H0=data(12,6,1)
def combined_energy(k,shift):
    return (2*H0(k)+(14+2*shift)*(v(1,k)**2+v(1,k-1)**2)
            -12*v(1,k)*v(1,k-1))
assert (combined_energy(7,0),combined_energy(8,0))==(46494,48746)
assert (combined_energy(7,1),combined_energy(8,1))==(47952,52516)
assert T(7)-T(11)>=abs(W(7,11))
print("combined-energy increases",46494,48746,47952,52516)

at,v,det,turb,T,W,H0=data(20,0,1)
assert F(5,6)*T(11)-T(12)==-2310
j,i=10,15
baseW=lambda j,i:v(0,j)*v(1,i)-v(1,j)*v(0,i)
K=lambda k:at(k)**2-at(k-2)*at(k+2)
extra=W(j,i)-baseW(j-1,i-1)-baseW(j+1,i+1)
assert (K(j)-K(i),extra)==(19404,34650)
assert (T(j)-T(i),W(j,i))==(108819,31500)
print("unchanged Sonin deficit",-2310,
      "signed insertion remainder",(19404,34650))

at,v,det,turb,T,W,H0=data(28,8,1)
j,k,i=14,16,24
aa=W(j,k)
zz=-W(j,k+1)
drop=T(k)-T(k+1)
assert aa>0 and zz>0 and drop>0
G=drop*(T(j)-T(i))-zz*(T(k)-T(i))-aa*(T(k+1)-T(i))
assert G==-740941856160
assert T(j)-T(i)-abs(W(j,i))==77554584
for sg in (-1,1):
    assert drop*(T(j)-T(i)+sg*W(j,i))==(
        zz*(T(k)-T(i)-sg*W(k,i))
        +aa*(T(k+1)-T(i)-sg*W(k+1,i))+G)
print("gluing correction",G)

at,v,det,turb,T,W,H0=data(32,8,2)
ii=(16,20,24,28)
j,k,l,m=ii
pf=W(j,k)*W(l,m)-W(j,l)*W(k,m)+W(j,m)*W(k,l)
assert pf==167057047945216
for u in range(4):
    for w in range(u+1,4):
        assert T(ii[u])-T(ii[w])>=abs(W(ii[u],ii[w]))
print("b=2 residual Pfaffian",pf,flush=True)

# Uniform band, including the sharper b=1 gap.
band=single=0
for N in (120,160,200,256,384,512):
    for de in sorted({0,2,4,2*(N//16),2*(N//8)}):
        sig=N+2
        if 4*de>sig:
            continue
        for B in range(1,7):
            at,v,det,turb,T,W,H0=data(N,de,B)
            gap=40 if B==1 else 36+12*B
            for j in range(N//2,N+1):
                xx=2*j-N-2*B
                if not (sig+8<=4*xx and 2*xx<=sig):
                    continue
                assert 135*sig*T(j)>=4*H0(j-B-1)
                for i in range(j+gap,N+B+1):
                    slack=T(j)-T(i)-abs(W(j,i))
                    assert 80*sig*slack>=H0(j-B-1)>0
                    assert N+B-j>=6
                    band+=1
                    single+=B==1
assert (band,single)==(234180,66500)
print("band intervals",band,"b=1",single,flush=True)

# An interior family for every b.
for B in range(1,13):
    a,e=72*B+72,56*B+56
    N=a+e
    de=a-e
    sig=N+2
    n,m=62*B+83,12*B+35
    j=(N+n-m)//2
    i=(N+n+m)//2+1
    at,v,det,turb,T,W,H0=data(N,de,B)
    assert n+m<N+2*B and N+B-j==40*B+40
    assert sig+8<=4*(n-m-2*B)
    assert 2*(n-m-2*B)<=sig and 4*de<=sig
    assert i-j==36+12*B
    assert 80*sig*(T(j)-T(i)-abs(W(j,i)))>=H0(j-B-1)>0
print("all-b interior family",12,flush=True)

# Independent SU(2) character multiplication.
def cg(i,j):
    return range(abs(i-j),i+j+1,2)

def mul(F0,G0):
    out=defaultdict(int)
    for (i,j),x in F0.items():
        for (k,l),y in G0.items():
            for u in cg(i,k):
                for w in cg(j,l):
                    out[u,w]+=x*y
    return {key:x for key,x in out.items() if x}

bridges=0
for a,e,B in ((3,3,1),(4,4,1),(5,5,1),(7,3,1),
              (4,2,2),(3,3,3),(9,3,1),(20,12,2)):
    N=a+e
    degree=N+2*B
    at,v,det,turb,T,W,H0=data(N,a-e,B)
    ff={(0,0):1}
    for factor,power in (({(1,0):1,(0,1):1},a),
                         ({(1,0):1,(0,1):-1},e),
                         ({(2,0):1,(0,2):1},B)):
        for _ in range(power):
            ff=mul(ff,factor)
    for n in range(3,degree+3):
        for m in range(3,n+1):
            if (N+n+m)%2:
                continue
            j=(N+n-m)//2
            i=(N+n+m)//2+1
            actual=(sum(ff.get((z,0),0) for z in cg(n,m)),
                    ff.get((n,m),0))
            assert actual==(T(j)-T(i),W(j,i))
            if B==1 and a==e:
                assert actual[0]>=abs(actual[1])
            bridges+=1
print("independent character bridges",bridges,flush=True)
print("PASS",flush=True)