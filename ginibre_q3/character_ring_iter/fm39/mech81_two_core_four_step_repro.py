import argparse
from math import comb
from fractions import Fraction as F
from collections import defaultdict
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

argparse.ArgumentParser(
    description="FM-MECH81 exact certificates and consumer checks."
).parse_args()

R,n,k,d=ring("n,k,d",QQ)
j=n+k
den={-2:(n+1)*(n+2),-1:n+1,0:R.one,1:R.one}
num={
    -2:(d*d-j*(n+1),-d*(j+1)),
    -1:(d,-j-1),0:(R.one,R.zero),1:(R.zero,R.one)
}
den[2]=j+2
num[2]=(-n,d)
for h in range(2,6):
    den[h+1]=den[h]*(j+h+1)
    num[h+1]=tuple(
        d*num[h][z]-(n-h+1)*(j+h)*num[h-1][z]
        for z in (0,1)
    )
H=(n+1)**2*(n+2)*(j+6)
for h in range(2,6):
    H*=(j+h)**2

def dot(h,l,weight=1):
    quot,rem=divmod(H,den[h]*den[l])
    assert not rem
    p,q=num[h]
    r,t=num[l]
    return [weight*quot*p*r,
            weight*quot*(p*t+q*r),
            weight*quot*q*t]

def plus(A,B):
    return [x+y for x,y in zip(A,B)]

rr=[R.zero]*3
ww=[R.zero]*3
terms=((0,0,1),(-1,-1,1),(1,1,1),
       (-2,0,-1),(0,2,-1),(-2,2,-1))
for h,sg in ((0,1),(4,-1)):
    for x,y,weight in terms:
        rr=plus(rr,dot(h+x,h+y,sg*weight))
for x in (-2,2):
    for y in (3,5):
        ww=plus(ww,dot(x,y,-1))
for x in (-1,1):
    for y in (2,6):
        ww=plus(ww,dot(x,y,1))
print("four-step forms built",flush=True)

P,u,v,w=ring("u,v,w",QQ)
def ev(poly,images):
    powers=[
        [P.one]+[x**h for h in range(1,poly.degree(z)+1)]
        for z,x in zip((n,k,d),images)
    ]
    out=P.zero
    for mon,cc in poly.items():
        term=P(cc)
        for po,h in zip(powers,mon):
            term*=po[h]
        out+=term
    return out

def pos(poly):
    assert poly and all(c>0 for c in poly.values())
    return len(poly)

def coefficient(poly,h):
    return R.from_dict({
        (i,j0,0):c for (i,j0,l),c in poly.items() if l==h
    })

def bern(parts):
    degree=len(parts)-1
    return [
        sum((QQ(comb(H0,h),comb(degree,h))*parts[h]
             for h in range(H0+1)),P.zero)
        for H0 in range(degree+1)
    ]

ABC={sg:plus(rr,[sg*z for z in ww]) for sg in (-1,1)}
FS={}
for sg,(A,B,C) in ABC.items():
    common=(n+1)**2*(d+sg*(2*n+k+2))**2
    for h in range(2,6):
        common*=(j+h)**2
    ff,rem=divmod(4*A*C-B*B,common)
    assert not rem and ff.degree(d)==10
    FS[sg]=ff

# Global positivity of the first principal coefficient.
for sg,(A,B,C) in ABC.items():
    nn,kk=u+5,v+2
    NN=2*nn+kk
    parts=[
        ev(coefficient(A,h),(nn,kk,P.zero))*(NN-4)**h
        for h in range(A.degree(d)+1)
    ]
    vals=bern(parts)
    assert all(q.get((0,0,0),0)>0 for q in vals)
    print("global A",sg,sum(pos(q) for q in vals),flush=True)

# Low-distance matrix region: two Bernstein expansions.
FE=(FS[-1]+FS[1])/2
FG=FS[-1]*FS[1]
kk=v+2
nn=5+(2*(kk+1)**2-5)*u
O=(nn-1)*(nn+kk+2)
for name,poly in (("even part",FE),("product",FG)):
    assert all(mon[2]%2==0 for mon in poly)
    degree=poly.degree(d)//2
    parts=[
        ev(coefficient(poly,2*h),(nn,kk,P.zero))*(4*O)**h
        for h in range(degree+1)
    ]
    total=0
    for pp in bern(parts):
        du=pp.degree(u)
        for H0 in range(du+1):
            qq=P.zero
            for (h,j0,l),cc in pp.items():
                if h<=H0:
                    qq+=P.from_dict({
                        (0,j0,l):cc*QQ(comb(H0,h),comb(du,h))
                    })
            total+=pos(qq)
    print("low matrix",name,total,flush=True)

# High-distance matrix region.
for sg in (-1,1):
    kk=v+2
    nn=2*(kk+1)**2+u
    LL=3*(kk+1)**2
    dc=(2*nn+kk)*(LL-2)
    parts=[
        ev(coefficient(FS[sg],h),(nn,kk,P.zero))
        *dc**h*LL**(10-h)
        for h in range(11)
    ]
    print("high matrix",sg,sum(pos(q) for q in bern(parts)),
          flush=True)

# Central ratio certificates.
for sg,abc in ABC.items():
    for odd in (False,True):
        ee=u+(3 if odd else 2)
        kk=v+2
        NN=3*ee*(kk+1)**2+w
        nn=(NN-kk)/2
        jj=(NN+kk)/2
        aa,bb,cc=[ev(x,(nn,kk,NN-2*ee)) for x in abc]
        if odd:
            dd=3*NN*kk
            lo=(kk+2)*(3*NN-(4*ee+2)*(kk+1))
            hi=3*NN*(kk+2)
        else:
            dd0=2*(NN-ee+2)-ee*kk**2
            dd=(jj+1)*dd0
            lo=nn*(dd0-4*ee*(kk+1))
            hi=nn*dd0
        vals=(
            aa*dd**2+bb*lo*dd+cc*lo**2,
            aa*dd**2+bb*(lo+hi)*dd/2+cc*lo*hi,
            aa*dd**2+bb*hi*dd+cc*hi**2
        )
        print("central",sg,odd,[pos(q) for q in vals],
              flush=True)

# Outer ratio certificates, with positive denominators cleared.
ee=u+2
tau=2*ee+1
LL=4*tau
xx=v
yy=v+w
nn=LL+xx*yy
kk=(2*ee-2)*LL+2*tau*(xx+yy)
dd=2*xx*yy+2*tau*(xx+yy)
degree=max(sum(mon) for pp in rr+ww for mon in pp)
assert degree==12
powers=[
    [P.one]+[pp**h for h in range(1,degree+1)]
    for pp in (nn,kk,dd,LL)
]
def lift(pp):
    ans=P.zero
    for mon,cc in pp.items():
        term=P(cc)
        for pow0,h in zip(powers,mon):
            term*=pow0[h]
        ans+=term*powers[3][degree-sum(mon)]
    return ans

R4,U,V,W,Z=ring("U,V,W,Z",QQ)
def outer_remainder(pp):
    reps=[R4.one,V]
    for h in range(2,pp.degree(v)+1):
        reps.append(
            -W*reps[-1]+(16*(2*U+5)+Z)*reps[-2]
        )
    return sum((
        cc*U**i*W**k0*reps[h]
        for (i,h,k0),cc in pp.items()
    ),R4.zero)

for sg,abc in ABC.items():
    aa,bb,cc=[lift(pp) for pp in abc]
    rn=nn
    rd=yy*(2*tau+xx)
    vals=(
        aa*rd**2,
        aa*rd**2+bb*rn*rd/2,
        aa*rd**2+bb*rn*rd+cc*rn**2
    )
    counts=[]
    for h,pp in enumerate(vals):
        if sg==-1 and h==1:
            pp=outer_remainder(pp)
        counts.append(pos(pp))
    print("outer",sg,counts,flush=True)

# Folded levels zero and one.
for sg,abc in ABC.items():
    for ee in (0,1):
        nn,kk=u+5,v+2
        NN=2*nn+kk
        jj=nn+kk
        aa,bb,cc=[ev(p0,(nn,kk,NN-2*ee)) for p0 in abc]
        p0=P.one if ee==0 else kk
        q0=nn if ee==0 else nn*(kk+2)
        val=(aa*((jj+1)*p0)**2
             +bb*((jj+1)*p0)*q0+cc*q0**2)
        if sg==-1 and ee==1:
            val,rem=divmod(val,(v+6)*NN*(NN+1))
            assert not rem
            val-=60*u**8*(v*v-2*u)**2
            val-=165*u**7*(2*v*v-11*u)**2
        print("strip",sg,ee,pos(val),flush=True)
print("ALL UNIFORM POLYNOMIAL CERTIFICATES PASS",flush=True)

def row(a,e):
    N=a+e
    out=[1]
    for h in range(N):
        z0,rem=divmod(
            (a-e)*out[-1]-(N-h+1)*(out[-2] if h else 0),
            h+1
        )
        assert not rem
        out.append(z0)
    return out

def data(a,e,b):
    N=a+e
    co=row(a,e)
    c=lambda h:co[h] if 0<=h<=N else 0
    def vh(h,k0):
        return sum(comb(h,l)*c(k0+h-2*l)
                   for l in range(h+1))
    T=lambda k0:sum(
        comb(b,h)*2**(b-h)*
        (vh(h,k0)**2-vh(h,k0-1)*vh(h,k0+1))
        for h in range(b+1)
    )
    W=lambda j0,i0:sum(
        comb(b,h)*2**(b-h)*
        (vh(h,j0)*vh(h+1,i0)-vh(h+1,j0)*vh(h,i0))
        for h in range(b+1)
    )
    return c,T,W

checks=0
minimum=None
coverage=[0,0,0,0]
for N0 in range(12,181):
    for e0 in range(2,(N0+1)//2):
        a0=N0-e0
        co=row(a0,e0)
        c=lambda h:co[h] if 0<=h<=N0 else 0
        BB=[c(h-1)+c(h+1) for h in range(N0+4)]
        AA=[-c(h-2)-c(h+2) for h in range(N0+4)]
        TT=[
            c(h)**2+c(h-1)**2+c(h+1)**2
            -c(h)*(c(h-2)+c(h+2))-c(h-2)*c(h+2)
            for h in range(N0+3)
        ]
        for j0 in range((N0+3)//2,N0-4):
            X0=2*j0-N0
            if X0<2:
                continue
            r0=N0-j0
            d0=a0-e0
            L0=3*(X0+1)**2
            if N0>=e0*L0:
                branch=0
            elif r0>=2*(X0+1)**2:
                branch=1
            elif d0*d0<=4*(r0-1)*(j0+2):
                branch=2
            else:
                branch=3
            if branch==0:
                rho=F(c(j0+1),c(j0))
                if e0%2:
                    hi=F(X0+2,X0)
                    lo=hi*(1-F((4*e0+2)*(X0+1),3*N0))
                else:
                    dd0=2*(N0-e0+2)-e0*X0*X0
                    hi=F(r0,j0+1)
                    lo=hi*(1-F(4*e0*(X0+1),dd0))
                assert lo<=rho<=hi
            if branch==3:
                p0,q0=abs(c(j0)),abs(c(j0+1))
                ss=d0*d0-4*(r0-1)*(j0+2)
                gap0=2*r0*p0-d0*q0
                assert c(j0)*c(j0+1)>0 and gap0>=0
                assert gap0*gap0>=ss*q0*q0
            coverage[branch]+=1
            RR=TT[j0]-TT[j0+4]
            WW=AA[j0]*BB[j0+4]-BB[j0]*AA[j0+4]
            slack=RR-abs(WW)
            assert slack>0
            rat=F(slack,TT[j0])
            if minimum is None or rat<minimum[0]:
                minimum=(rat,a0,e0,j0)
            checks+=1
assert checks==434690
print("target search",checks,"coverage",coverage,
      "minimum",minimum,flush=True)

def eval3(pp,n0,k0,d0):
    return sum(cc*n0**h*k0**l*d0**s0
               for (h,l,s0),cc in pp.items())

for rr0,xx0,dd0,value in (
    (5,19,25,F(-174394727488,5076990009)),
    (249,2,496,F(-14763544884246518929837321,
                  4966152000152551111684000000))
):
    aa,bb,cc=[
        eval3(pp,rr0,xx0,dd0)/eval3(H,rr0,xx0,dd0)
        for pp in ABC[1]
    ]
    assert 4*aa*cc-bb*bb==value<0
print("arbitrary-state obstructions PASS",flush=True)

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

def mul(f,g):
    out=defaultdict(int)
    for (a,b),x in f.items():
        for (c,d),y in g.items():
            for i in cg(a,c):
                for j0 in cg(b,d):
                    out[i,j0]+=x*y
    return {ij:x for ij,x in out.items() if x}

bridges=0
for a0,e0 in ((9,3),(10,4),(12,4),(11,5),(13,3),(14,6)):
    N0=a0+e0
    ff={(0,0):1}
    for fac,power in (
        ({(1,0):1,(0,1):1},a0),
        ({(1,0):1,(0,1):-1},e0),
        ({(2,0):1,(0,2):1},1)
    ):
        for _ in range(power):
            ff=mul(ff,fac)
    c,T,W=data(a0,e0,1)
    for j0 in range((N0+3)//2,N0-4):
        X0=2*j0-N0
        if X0<2:
            continue
        label=X0+3
        RR=sum(ff.get((h,0),0) for h in cg(label,3))
        WW=ff.get((label,3),0)
        assert (RR,WW)==(T(j0)-T(j0+4),W(j0,j0+4))
        hh=eval3(H,N0-j0,X0,a0-e0)
        for sg,abc in ABC.items():
            aa,bb,cc=[
                eval3(pp,N0-j0,X0,a0-e0) for pp in abc
            ]
            assert (
                aa*c(j0)**2+bb*c(j0)*c(j0+1)
                +cc*c(j0+1)**2==hh*(RR+sg*WW)
            )
        for eta in (-1,1):
            eps=(-1)**e0*eta
            word=mul(
                mul(ff,{(label,0):1,(0,label):eps}),
                {(3,0):1,(0,3):eta}
            )
            assert word.get((0,0),0)==2*(RR+eta*WW)>0
            bridges+=1
assert bridges==34
print("independent full-word bridges",bridges,flush=True)

for aa,ee,bb,j0,g0,expected in (
    (10,6,1,10,5,(10935,-210,10725)),
    (9,3,2,7,4,(8452,-788,7664))
):
    c,T,W=data(aa,ee,bb)
    RR=T(j0)-T(j0+g0)
    WW=W(j0,j0+g0)
    assert (RR,WW,RR-abs(WW))==expected

# Bounded extension search; not a uniform theorem.
count=0
best=None
for N0 in range(12,65):
    for e0 in range(2,(N0+1)//2):
        a0=N0-e0
        co=row(a0,e0)
        c=lambda h:co[h] if 0<=h<=N0 else 0
        for b0 in range(1,7):
            VV=[
                [sum(comb(h,l)*c(k0+h-2*l)
                     for l in range(h+1))
                 for k0 in range(-1,N0+b0+4)]
                for h in range(b0+2)
            ]
            vh=lambda h,k0:VV[h][k0+1]
            wt=[comb(b0,h)*2**(b0-h) for h in range(b0+1)]
            TT=[
                sum(wt[h]*(vh(h,k0)**2
                    -vh(h,k0-1)*vh(h,k0+1))
                    for h in range(b0+1))
                for k0 in range(N0+b0+2)
            ]
            for j0 in range((N0+1)//2,N0+b0-5):
                for gap in range(4,min(14,N0+b0-j0+1)):
                    i0=j0+gap
                    n0=2*j0-N0+gap-1
                    m0=gap-1
                    if n0<m0 or n0+m0>=N0+2*b0:
                        continue
                    RR=TT[j0]-TT[i0]
                    WW=sum(
                        wt[h]*(vh(h,j0)*vh(h+1,i0)
                               -vh(h+1,j0)*vh(h,i0))
                        for h in range(b0+1)
                    )
                    slack=RR-abs(WW)
                    assert slack>0
                    rat=F(slack,TT[j0])
                    count+=1
                    if best is None or rat<best[0]:
                        best=(rat,a0,e0,b0,j0,gap)
assert count==979270
print("bounded extension search",count,"minimum",best,
      flush=True)
print("PASS",flush=True)