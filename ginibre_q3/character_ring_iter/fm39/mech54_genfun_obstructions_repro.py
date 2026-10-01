import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb

argparse.ArgumentParser(
    description="FM-MECH54: generating identities and scoped obstructions"
).parse_args()

def C(n,k):
    return comb(n,k) if 0 <= k <= n else 0

def cat(j):
    return C(2*j,j)//(j+1)

@lru_cache(None)
def row(e,a):
    c=[1]
    for sign,count in ((1,a),(-1,e)):
        for _ in range(count):
            c=[(c[j] if j<len(c) else 0)
               +sign*(c[j-1] if j else 0)
               for j in range(len(c)+1)]
    return tuple(c)

@lru_cache(None)
def profile(e,a,b):
    T=e+a
    D=T+2*b
    mon=[0]*(D+1)
    c=row(e,a)
    for j in range(0,T+1,2):
        for l in range(b+1):
            for h in range(l+1):
                mon[T-j+2*(l-h)]+=(
                    c[j]*C(b,l)*(-2)**(b-l)*C(l,h)*cat(j//2+h)
                )
    return tuple(
        sum(mon[m]*(C(m,(m-p)//2)-C(m,(m-p)//2-1))
            for m in range(p,D+1,2))
        for p in range(D+1)
    )

def pmul(P,R):
    out={}
    for (i,j),v in P.items():
        for (k,l),w in R.items():
            key=(i+k,j+l)
            out[key]=out.get(key,0)+v*w
    return {key:v for key,v in out.items() if v}

GENERATORS={
    "-":{(0,0):1,(2,0):1,(1,1):-1},
    "+":{(0,0):1,(2,0):1,(1,1):1},
    "B":{(0,0):1,(4,0):1,(2,2):1},
    "V":{(0,2):1,(0,0):-1},
}

@lru_cache(None)
def power(name,n):
    out={(0,0):1}
    for _ in range(n):
        out=pmul(out,GENERATORS[name])
    return out

@lru_cache(None)
def gf(e,a,b,h=0):
    P=pmul(pmul(power("-",e),power("+",a)),power("B",b))
    P=pmul(P,power("V",2*h))
    H=[0]*(e+a+2*b+1)
    for (i,j),v in P.items():
        if j%2==0:
            assert i%2==0
            H[i//2]+=v*cat(j//2)
    return tuple(
        (H[d] if d<len(H) else 0)-(H[d-1] if d else 0)
        for d in range(len(H)+1)
    )

def zmul(P,R):
    out=[0]*(len(P)+len(R)-1)
    for i,v in enumerate(P):
        for j,w in enumerate(R):
            out[i+j]+=v*w
    return out

def zpower(P,n):
    out=[1]
    for _ in range(n):
        out=zmul(out,P)
    return out

bridges=0
for e in range(7):
    for a in range(7):
        for b in range(5):
            D=e+a+2*b
            G=gf(e,a,b)
            assert all(G[j]==-G[D+1-j] for j in range(D+2))
            for d in range(D//2+1):
                assert G[d]==profile(e,a,b)[D-2*d]
                bridges+=1
print("Generating/direct Catalan bridges:",bridges)

signed_checks=0
ratio_checks=0
for b in range(4):
    for e in range(b,b+5):
        for a in range(b,b+5):
            G=gf(e,a,b)
            rhs=[0]*len(G)
            for h in range(b+1):
                term=zmul(zpower([1,1,1],2*b-2*h),
                           gf(e-b,a-b,0,h))
                for j,v in enumerate(term):
                    rhs[j+2*h]+=(-1)**h*C(b,h)*v
            assert tuple(rhs)==G
            signed_checks+=1

            denominator=gf(e-b,a-b,0)
            coeff=[]
            for d in range(4):
                value=G[d] if d<len(G) else 0
                value-=sum(
                    coeff[h]*
                    (denominator[d-h] if d-h<len(denominator) else 0)
                    for h in range(d)
                )
                coeff.append(value)
            T=e+a
            u=(a-e)**2
            expected=[1,2*b,2*b*b,
                      Q(b*(3*T+8*b*b-6*b-3*u+4),6)]
            assert coeff==expected
            ratio_checks+=1
print("Signed background identities:",signed_checks)
print("Generating-ratio coefficient identities:",ratio_checks)

den=gf(0,8,0)
num=gf(3,11,3)
rat=[]
for d in range(4):
    rat.append(num[d]-sum(rat[h]*den[d-h] for h in range(d)))
assert rat==[1,6,18,-46]
print("Ratio witness (e,a,b)=(3,11,3):",rat)

minima=[]
for L in range(8,21,2):
    vals=[]
    for e in range(L//2+1):
        f=profile(e,L-e,0)
        vals.append(17*f[8]-19*sum(f[10::2]))
    assert min(vals)>=0
    minima.append((L,min(vals)))

target=profile(6,8,3)
assert [target[p] for p in range(8,21,2)]==[
    2120,1295,613,214,58,11,1
]
assert 17*target[8]-19*sum(target[10::2])==-5608
print("Whole-tail separator minima:",minima)
print("Whole-tail separator on g_(6,8,3):",-5608)

# Exact polynomials in the imbalance xi.
def padd(P,R,scale=1):
    out=list(P)+[Q(0)]*max(0,len(R)-len(P))
    for j,v in enumerate(R):
        out[j]+=scale*v
    return out

def craw(T,n):
    c=[[Q(1)],[Q(0),Q(1)]]
    for j in range(1,n):
        nxt=padd([Q(0)]+c[j],c[j-1],-(T-j+1))
        c.append([v/(j+1) for v in nxt])
    return c

@lru_cache(None)
def mu(j,l):
    return sum(C(l,h)*(-2)**(l-h)*cat(j+h) for h in range(l+1))

def fpoly(T,b,d):
    c=craw(T,2*d)
    out=[Q(0)]*(d+1)
    for j in range(d+1):
        factor=sum(
            C(b,l)*mu(j,l)*
            (C(T+2*b-2*j-2*l,d-j-l)
             -C(T+2*b-2*j-2*l,d-j-l-1))
            for l in range(min(b,d-j)+1)
        )
        for h in range(j+1):
            out[h]+=c[2*j][2*h]*factor
    return out

def peval(P,u):
    return sum(v*u**j for j,v in enumerate(P))

weights=[Q(1),Q(6),Q(54,7),Q(88,5),Q(813,5),Q(498),Q(-216)]
difference=fpoly(20,3,6)
for h,w in enumerate(weights):
    difference=padd(difference,fpoly(14,0,6-h),-w)
assert not any(difference)

for e in range(3,11):
    a=20-e
    actual=profile(e,a,3)[14]
    reduced=sum(
        weights[h]*profile(e-3,a-3,0)[2+2*h]
        for h in range(7)
    )
    assert actual==reduced==peval(fpoly(20,3,6),(a-e)**2)
print("Shifted-base polynomial identity:",list(map(str,weights)))
print("Shifted-base actual-row bridges:",8)

xis=[1,3,5,7]
values=[]
rows=[]
for xi in xis:
    ep=(7-xi)//2
    ap=(7+xi)//2
    c=row(ep,ap)
    actual=profile(ep+3,ap+3,3)[7]
    formula=(
        Q(4634,5)+102*c[1]**2
        +Q(106,5)*c[2]**2-Q(2,5)*c[3]**2
    )
    assert actual==formula==peval(fpoly(13,3,6),xi*xi)
    values.append(actual)
    rows.append(c)

assert values==[1216,1856,5184,14784]
dual=[-5,9,-5,1]
square_duals=[]
for j in range(8):
    num=sum(dual[i]*rows[i][j]**2 for i in range(4))
    assert num==1280*(j in (3,4))
    square_duals.append(num//1280)
assert Q(sum(dual[i]*values[i] for i in range(4)),1280)==Q(-2,5)
print("Diagonal-square actual-row values:",values)
print("Dual values on the eight row squares:",square_duals)
print("Dual value on the target:",Q(-2,5))

vals=[profile(7,7,b)[2+2*b] for b in range(3,10)]
newton=[]
while vals:
    newton.append(vals[0])
    vals=[vals[j+1]-vals[j] for j in range(len(vals)-1)]
assert newton==[1832,798,404,274,-32,88,132]
print("Residual Newton coefficients:",newton)
print("PASS")
