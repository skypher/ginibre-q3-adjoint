
from math import comb, sqrt, isqrt
from collections import Counter

SCALE = 10**24

def row(a,e):
    N=a+e
    return [sum((-1)**u*comb(e,u)*comb(a,k-u)
                for u in range(e+1) if 0<=k-u<=a)
            for k in range(N+1)]

def at(c,k):
    return c[k] if 0<=k<len(c) else 0

def ceil_sqrt_int(x):
    if x<=0: return 0
    y=x*SCALE*SCALE
    r=isqrt(y)
    return r if r*r==y else r+1

def ceil_sqrt_ratio(num,den):
    if num<=0: return 0
    y=num*SCALE*SCALE
    r=isqrt(y//den)
    return r if r*r*den==y else r+1

def metric_chord_upper(x,y,D,N,V):
    C=y-x
    product=D(x)*D(y)
    ld=ceil_sqrt_int((C+1)**2*product)
    h=2*x+C-N
    Rp=C+h
    Rm=abs(C-h)
    K=(N+2-C)**2-h*h
    if K<0 or Rp*Rp>=4*V:
        return ld
    U=K*product
    A=4*V+Rp*Rm
    B=2*(Rp+Rm)
    Q=A*A-B*B*V
    assert Q>0
    sqrtV_upper=ceil_sqrt_int(V)
    # The squared metric bound is U/(A-B*sqrt(V)).
    # Rationalize it, then use upper bounds for both square roots.
    metric=ceil_sqrt_ratio(U*(A*SCALE+B*sqrtV_upper),Q*SCALE)
    return min(ld,metric)

def exact_total_upper(s,D,N,V,j,i):
    x=D(j)*D(i-1)
    y=D(j+1)*D(i)
    ld=s*(ceil_sqrt_int(x)+ceil_sqrt_int(y))
    metric=(metric_chord_upper(j,i-1,D,N,V)
            +metric_chord_upper(j+1,i,D,N,V))
    return min(ld,metric)

def float_chord(x,y,D,N,V):
    C=y-x
    ld=(C+1)*sqrt(max(D(x)*D(y),0))
    h=2*x+C-N
    Rp=C+h
    Rm=abs(C-h)
    K=(N+2-C)**2-h*h
    two_sqrtV=2*sqrt(V)
    if Rp<two_sqrtV and K>=0:
        m=sqrt(K*D(x)*D(y)/
               ((two_sqrtV-Rp)*(two_sqrtV-Rm)))
        return min(ld,m)
    return ld

def scan(amin,amax,scaled,tol,arange=None):
    stats=Counter()
    uncertified=[]
    for a in (arange if arange is not None else range(amin,amax+1)):
        for e in range(amin,amax+1):
            N=a+e
            V=(a+1)*(e+1)
            c=row(a,e)
            C=lambda k:at(c,k)
            Dint=lambda k:C(k)**2-C(k-1)*C(k+1)
            scale=max(Dint(k) for k in range(N+2)) if scaled else 1
            Dfloat=lambda k:Dint(k)/scale

            for j in range((N+1)//2,N+1):
                for i in range(j+1,N+2):
                    q=i-j-1
                    s=i-j
                    stats["pairs"]+=1
                    if abs(a-e)<=1 or q<=1 or (j,i)==(N-1,N+1):
                        stats["skipped"]+=1
                        continue
                    if not scaled and min(a,e)<=2:
                        stats["skipped"]+=1
                        continue

                    Lfloat=Dfloat(j)-Dfloat(i)
                    ld=s*(sqrt(Dfloat(j)*Dfloat(i-1))
                          +sqrt(Dfloat(j+1)*Dfloat(i)))
                    b=min(ld,
                          float_chord(j,i-1,Dfloat,N,V)
                          +float_chord(j+1,i,Dfloat,N,V))
                    if Lfloat>=b*(1+tol):
                        stats["float accepted"]+=1
                        upper=exact_total_upper(s,Dint,N,V,j,i)
                        if (Dint(j)-Dint(i))*SCALE < upper:
                            uncertified.append((a,e,j,i))
                        else:
                            stats["integer certified"]+=1
    return stats,uncertified


import time
t0=time.time(); last=t0; T=Counter(); B=[]
for a0 in range(3,81):
    st_,bad_=scan(3,80,True,1e-9,arange=[a0]); T.update(st_); B+=bad_
    if time.time()-last>30 or a0==80:
        last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat a=%d/80'%a0,dict(T),'uncertified',len(B),flush=True)
print('done',dict(T),'uncertified',len(B),B[:5],'elapsed',round(time.time()-t0,1))
