"""FM-MECH74: larger-anchor theorem and boundary phase criterion.
Run from the repository root with python3 -u; no files are written.
"""
import argparse
from fractions import Fraction as F
from math import comb, factorial, isqrt
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n",type=int,default=160,
                help="bounded diagnostic census; not part of the proof")
args=ap.parse_args()

# 1. The exact binomial-energy transfer.
si,de,X,u,v=S.symbols("si de X u v", real=True)
C=si*si-de*de
V=C-X*X+2*(si-X)
up=(de*u-V*v)/(si+X+2)
vp=(u+de*v)/(si+X+2)
assert S.cancel(u-(de*up+V*vp)/(si-X))==0
assert S.cancel(v-(-up+de*vp)/(si-X))==0
En=u*u+V*v*v
Vn=C-(X+2)**2+2*(si-X-2)
assert S.cancel(up*up+Vn*vp*vp
                -(si-X)*En/(si+X+2)+4*(X+2)*vp*vp)==0
assert S.expand(de*de+V-(si-X)*(si+X+2))==0

# The actual change from (c_k,B_k), and its recurrence.
p,q=S.symbols("p q")
pn=(2*de*p-(si-X)*q)/(si+X)
pnn=(2*de*pn-(si-X-2)*p)/(si+X+2)
bb=q+pn
uu=si*p-de*bb/2
vv=(si*bb/2-de*p)/X
uun=si*pn-de*(p+pnn)/2
vvn=(si*(p+pnn)/2-de*pn)/(X+2)
assert S.cancel(uun-(de*uu-V*vv)/(si+X+2))==0
assert S.cancel(vvn-(uu+de*vv)/(si+X+2))==0
print("TRANSFER IDENTITIES: PASS")

# 2. Maximal turn for [[d,-b],[c,d]], 0 <= b <= c.
d,b,c=S.symbols("d b c",positive=True)
aa=(b+c)/(2*d); bb=(c-b)/(2*d)
assert S.cancel((1+aa*aa+aa*bb)**2-(1+aa*aa-bb*bb)
                -c*c*((b+c)**2+4*d*d)/(4*d**4))==0
rho=S.symbols("rho")
bound=(F(151,100)**2*(1-rho*rho)
       -F(65,64)**2*F(16,25)
       -F(65,64)*F(1,15)*(1-F(3,5)*rho))
assert S.Poly(bound,rho).coeff_monomial(rho*rho)<0
assert bound.subs(rho,0)>0
assert bound.subs(rho,F(5,6))==F(3869,1440000)
assert F(101,1000)-F(101,1000)**3/6>F(1,10)
assert F(151,50)+F(101,1000)<F(157,50)
assert 16*(F(1,5)-F(1,375))-F(4,239)>F(157,50)
print("LONG INNER ARC: kappa*(y-x)>2, CONSTANTS PASS")

# 3. The (BE) certificate at y=x+2/kappa.
x,t,U,T=S.symbols("x t U T")
y=x+2*t
q=4*x+4*t
poly=S.expand((1-y)*(1-x*x)*sum(q**h/S.factorial(h) for h in range(9))
              -(1+y)*(1-x*x+t))
pp=S.Poly(S.expand(poly.subs({
    x:S.Rational(3,5)+U/S.Integer(10),t:T/S.Integer(15)})),U,T)
du,dt=pp.degree(U),pp.degree(T)
terms=dict(pp.terms())
bern=[]
for i in range(du+1):
    for j in range(dt+1):
        val=sum(z*S.Rational(comb(i,h),comb(du,h))*
                    S.Rational(comb(j,k),comb(dt,k))
                for (h,k),z in terms.items() if h<=i and k<=j)
        bern.append(val)
assert (du,dt,len(bern))==(11,9,120)
assert min(bern)==S.Rational(61336546830917,80731054687500)
print("BE BERNSTEIN",du,dt,len(bern),"MIN",min(bern))

# 4. Checkpoint z in [93/100,769/800].
def exp_lower(q,n=24):
    return sum(q**h/F(factorial(h)) for h in range(n+1))
for z in (F(93,100),F(769,800)):
    exponent=15*(z*z-F(7,10)**2)
    prefactor=(1+z)/((1+F(7,10))*(1-F(7,10))*(1-z))
    assert exp_lower(exponent)>5*prefactor
assert 4*F(93,100)/(1-F(93,100)**2)**2>30
assert F(1,5)*F(2,5)*F(7,100)<F(1,100)
assert F(1,5)*F(16,15)<1
margin=(1-F(1,100))**2-4*F(33,32)**2/F(5)
assert margin==F(20691,160000)
print("CHECKPOINT MARGIN",margin)

def data(a,e):
    N=a+e;si=N+2;de=a-e;C=4*(a+1)*(e+1)
    row=[1]
    for k in range(N):
        val,rem=divmod(de*row[k]-(N-k+1)*(row[k-1] if k else 0),k+1)
        assert rem==0
        row.append(val)
    c=lambda k:row[k] if 0<=k<=N else 0
    B=[c(k-1)+c(k+1) for k in range(N+2)]
    D=[c(k)**2-c(k-1)*c(k+1) for k in range(N+2)]
    H=[si*(c(k)**2+c(k-1)**2)-2*de*c(k)*c(k-1)
       for k in range(N+2)]
    beta=[comb(si,k+1) for k in range(N+2)]
    return N,si,de,C,c,B,D,H,beta

def be(N,si,C,j,i,beta):
    X,Y=2*j-N,2*i-N
    P=C-X*X;V=P+2*(si-X)
    aa=P*beta[j]-V*beta[i]
    bb=Y*(P*beta[j]+V*beta[i])
    return aa>=0 and C*aa*aa>=bb*bb

def joint(N,si,C,j,i,D,H):
    Y=2*i-N
    return C*(si+Y)**2*(D[j]-D[i])**2-4*Y*Y*H[j]*H[i]

# 5. Bounded checks of assertions used in the new anchor theorem.
rows=anchors=pairs=long_inner=early=tail=0
samples={(N-e,e) for N in range(80,args.max_n+1)
         for e in range(7,(N-2)//2+1)}
samples.update((a,e) for a in (400,1000)
               for e in (7,8,12,20,40,60))
for a,e in sorted(samples):
    N=a+e;si=N+2;C=4*(a+1)*(e+1)
    if C<4096 or C<30*(si+1):continue
    js=[j for j in range(N//2+1,N-3)
        if 25*(2*j-N)**2>9*C and 100*(2*j-N)**2<=49*C
        and 2*(2*j-N)<=si]
    if not js:continue
    rows+=1
    N,si,de,C,c,B,D,H,beta=data(a,e)
    m=(N+isqrt(C))//2
    while 10000*(2*m-N)**2>=8649*C:m-=1
    m+=1
    for j in js:
        anchors+=1;X=2*j-N
        assert m>j and (2*m-N)**2<C
        seen=False
        for i in range(j+1,N+2):
            Y=2*i-N
            W=c(j)*B[i]-B[j]*c(i)
            seen |= W<0 or (W==0 and c(j)*c(i)+B[j]*B[i]<0)
            if i-j<4:continue
            assert D[j]-D[i]>=abs(W),(a,e,j,i)
            pairs+=1
            if seen and Y*Y<C:
                assert C*(Y-X)**2>16*(si+1)**2
                long_inner+=1
            if seen and i<m:
                assert be(N,si,C,j,i,beta),(a,e,j,i)
                early+=1
            if i>=m:
                assert joint(N,si,C,j,i,D,H)>=0,(a,e,j,i)
                tail+=1
print("ANCHOR THEOREM",rows,anchors,pairs,
      "INNER LONG",long_inner,"EARLY BE",early,"TAIL J",tail)

# 6. Boundary direction and the parameter-only no-long-arc region.
nr=na=nn=nl=no=0
for a in (80,100,160,240,400,1000):
    for e in (7,8,9,12,15,20,30,40,60,90):
        if a<e+2:continue
        N,si,de,C,c,B,D,H,beta=data(a,e)
        if C<4096 or C<30*(si+1):continue
        nr+=1
        K=(N+isqrt(C))//2
        while (2*K-N)**2<C:K+=1
        sg=(-1)**e
        for k in range(K,N+2):
            X=2*k-N
            uu=sg*(2*si*c(k)-de*B[k])
            vv=sg*(si*B[k]-2*de*c(k))
            assert uu<=0 and vv>0
            assert X*X*uu*uu>=(X*X-C)*vv*vv
            no+=1
        for j in range(N//2+1,K):
            X=2*j-N
            if 25*X*X<=9*C or 2*X>si:continue
            na+=1;s=K-j
            V=C-X*X+2*(si-X)
            SV=s*V-2*(X+1)*s*(s-1)-2*s*(s-1)*(2*s-1)//3
            assert SV==sum(C-(X+2*h)**2+2*(si-X-2*h) for h in range(s))
            small=10000*s*SV<=24649*de*de
            seen=False
            for i in range(j+1,N+2):
                W=c(j)*B[i]-B[j]*c(i)
                seen |= W<0 or (W==0 and c(j)*c(i)+B[j]*B[i]<0)
            if small:
                assert not seen,(a,e,j,K)
                nn+=1
            if seen:
                assert not small
                nl+=1
assert (nr,na,nn,nl,no)==(53,1162,265,698,5731)
print("BOUNDARY PHASE",nr,na,"NO-LONG",nn,"LONG",nl,"OUTER",no)

# The earlier short-arc turning control now has a uniform certificate.
a,e,j=160,9,120
N,si,de,C,c,B,D,H,beta=data(a,e)
K=125;X=2*j-N;s=K-j
V=C-X*X+2*(si-X)
SV=s*V-2*(X+1)*s*(s-1)-2*s*(s-1)*(2*s-1)//3
assert (C,X,s,SV,de*de)==(6440,71,5,4995,22801)
assert 10000*s*SV<24649*de*de
assert joint(N,si,C,j,K,D,H)<0
assert all(c(j)*B[i]-B[j]*c(i)>=0 for i in range(j+1,N+2))
print("SHORT TURNING CONTROL: PARAMETER CRITERION PASS; J(K)<0")

# Outside the two new uniform sufficient regions.
a,e,j,i=78,12,68,74
N,si,de,C,c,B,D,H,beta=data(a,e)
X,Y=2*j-N,2*i-N
K=(N+isqrt(C))//2
while (2*K-N)**2<C:K+=1
s=K-j;V=C-X*X+2*(si-X)
SV=s*V-2*(X+1)*s*(s-1)-2*s*(s-1)*(2*s-1)//3
assert (N,si,C,X,Y,K)==(90,92,4108,46,58,78)
assert 100*X*X>49*C and 4*s*SV>10*de*de
assert X*X<4*(e-1)*(a+4)
assert c(j)*B[i]-B[j]*c(i)<0
slack=D[j]-D[i]-abs(c(j)*B[i]-B[j]*c(i))
assert slack==21039816540244478436889833766635
assert joint(N,si,C,j,i,D,H)>0
print("RESIDUAL SAMPLE",(a,e,j,i),"POSITIVE SLACK",slack)
print("PASS")

