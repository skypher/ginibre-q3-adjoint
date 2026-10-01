from fractions import Fraction as Q
from math import comb,factorial,isqrt
import sympy as S

# Critical-sector quartic certificate.
v,w=S.symbols("v w")
A=8/(4-v)+4/(2+v)
C=(S.Rational(2,3)/(4-v)+S.Rational(4,3)/(2+v)
   +(7+2*v)/((4-v)*(S.Rational(49,16)-v))
   +8/((2+v)*(1+v)))
num,den=S.fraction(S.factor(S.Rational(3,8)*A*A-C))
num,den=-num,-den
P=S.Poly(S.expand(num.subs(v,S.Rational(9,4)*w)),w)
N=P.degree()
bc=[sum(P.nth(j)*S.Rational(comb(i,j),comb(N,j))
        for j in range(i+1)) for i in range(N+1)]
assert min(bc)>0
assert Q(12,11)**7<Q(11,8)**2
assert Q(12,11)**5<Q(5,4)**2
assert Q(12,11)**3<Q(8,7)**2
assert Q(1,12)+Q(495,256)<Q(33,16)
assert Q(64,35)<Q(15,8)
assert Q(9,5)*Q(53,45)*15<32

UNIT=2**100
def root(x):
    j=isqrt(x.numerator*UNIT*UNIT//x.denominator)
    return Q(j,UNIT),Q(j+1,UNIT)
def exp_lo(x,n=24):
    term=out=Q(1)
    for j in range(1,n+1):
        term*=x/Q(j)
        out+=term
    return out
def exp_hi(x,n=24):
    assert x<n+2
    return exp_lo(x,n)+x**(n+1)/factorial(n+1)/(1-x/Q(n+2))

def critical_bound(b,M):
    r=root(Q(b))[1]
    tau=Q(2*M*M,31*b)
    return exp_hi(tau,48)*(
        Q(33,16*b)+Q(15,8*M)
        +32*b*r/exp_lo(Q(31*b,36),48)
        +(2*r+Q(1,isqrt(b)))/exp_lo(Q(31*b,32),48)
        +200*b*b*Q(5,12)**b)

for b,M,limit in ((16,9,530),(16,17,910),(32,31,868),
                  (64,33,267),(64,49,794),(128,65,379)):
    assert critical_bound(b,M)<Q(limit,1000)
tail=(Q(32,12)*27**2*Q(21,4)/exp_lo(Q(93,4),48)
      +Q(963,40)/exp_lo(Q(837,32),48)
      +Q(50,3)*27**3*Q(5,12)**27)
assert tail<Q(1,64)
print("Critical-sector bounds: PASS")

# Maximizer identities.
y,a,e,b=S.symbols("y a e b")
T=a+e
Lp=-e/(2-y)+a/(2+y)+2*b*y/(2+y*y)
Ay=e/(2-y)+a/(2+y)+4*b/(2+y*y)
P=(T+2*b)*y**3-2*(a-e)*y*y+(2*T-8*b)*y-4*(a-e)
assert S.factor(P+(4-y*y)*(2+y*y)*Lp)==0
assert S.factor(Ay-(T/2+b)-2*b/(2+y*y)+y*Lp/2)==0
assert S.factor((2-y)*Ay-2*e-2*b*(2-y)**2/(2+y*y)
                -(2-y)*Lp)==0

# Exact no-margin witnesses.
def choose(n,j):
    return comb(n,j) if 0<=j<=n else 0
for m,digits in ((20,9),(24,21),(32,57)):
    r=m*m+m-3
    M=2*m*m+1
    D=2*r+6
    F=0
    for j in range(m+1):
        for l in range(min(3,m-j)+1):
            mu=sum(choose(l,h)*(-2)**(l-h)
                   *choose(2*(j+h),j+h)//(j+h+1)
                   for h in range(l+1))
            n=D-2*j-2*l
            q=m-j-l
            F+=((-1)**j*choose(r,j)*choose(3,l)*mu
                *(choose(n,q)-choose(n,q-1)))
    bound=(Q(4*(r+7)**2,3*M)*choose(2*r+13,m)*Q(3,4)**r)
    assert F>0 and bound<Q(1,10**digits)
assert Q(64**8,factorial(8)*864*512**2)>1
print("No-margin witnesses and escaping family: PASS")

# Constants used in the endpoint interval bounds.
def atan_iv(d):
    q=sum((Q((-1)**j,(2*j+1)*d**(2*j+1))
           for j in range(24)),Q(0))
    return q,q+Q(1,49*d**49)
p5,p239=atan_iv(5),atan_iv(239)
pi_lo=16*p5[0]-4*p239[1]
pi_hi=16*p5[1]-4*p239[0]
assert Q(177,100)**2<pi_lo<pi_hi<Q(9,5)**2
assert Q(5,48)+Q(165,32)*Q(11,18)+Q(8,7)<5
assert Q(5,2)**3<16

def P_hi(x):
    lo,hi=root(x)
    return (2*hi+1/lo)*Q(100,177)/exp_lo(x)
def beta_x_int(x,n):
    q=Q(factorial(n-1))
    for j in range(n):q/=x+j
    return q
def negative(e,a,b):
    first=beta_x_int(Q(a+e,2)+1,b+1)
    if a%2==e%2==0:
        A,E=a//2,e//2
        q=Q(factorial(2*A)*factorial(2*E),
            4**(A+E)*factorial(A)*factorial(E)*factorial(A+E))
        return Q(36,5)*first*q
    if a%2:
        second=beta_x_int(Q(e+1,2),(a+1)//2)
    else:
        second=beta_x_int(Q(a+1,2),(e+1)//2)
    return Q(400,177)*first*second

def certify(e,a,b,p,N,limit):
    M=p+1
    delta=Q(1,512)
    dy=(4-2*delta)/N
    Gticks=Eticks=0
    for j in range(N):
        lo=-2+delta+j*dy
        hi=lo+dy
        al=Q(0) if lo<=0<=hi else min(abs(lo),abs(hi))
        ah=max(abs(lo),abs(hi))
        vl,vh=al*al,ah*ah
        sl,sh=4-vh,4-vl
        d1l,d1h=2-hi,2-lo
        d2l,d2h=2+lo,2+hi
        Al=Q(e)/d1h+Q(a)/d2h+Q(4*b)/(2+vh)
        Ah=Q(e)/d1l+Q(a)/d2l+Q(4*b)/(2+vl)
        Llo=min(Q(1,4),(2-ah)/2)
        Lhi=min(Q(1,4),(2-al)/2)
        c=1-Lhi/3
        k=max(Q(2,5),1-Q(5,12)*(2-al))
        C=(Q(e)/d1l+Q(a)/d2l)/12+Q(4*b,3)/(2+vl)
        C+=Q(e)/d1l**2+Q(a)/d2l**2+Q(16*b)/(2+vl)**2
        J=min(Q(M),2/root(sl)[0])
        er=Q(1,4)/Al*root(1/c**5)[1]
        er+=Q(15,4)*C/Al**2*root(1/c**7)[1]
        er+=J/M*root(1/c**3)[1]
        er+=(1+J/M)*root(1/k**3)[1]*P_hi(k*Al*Llo)
        er+=P_hi(Al*Llo)
        hl=(1-hi/2)**e*(1+lo/2)**a*(1+vl/2)**b
        hh=(1-lo/2)**e*(1+hi/2)**a*(1+vh/2)**b
        gl=dy*hl*root(sl/Ah**3)[0]/exp_hi(Q(M*M)/(4*Al))
        eh=dy*hh*root(sh/Al**3)[1]*er
        Gticks+=gl.numerator*UNIT//gl.denominator
        Eticks+=(eh.numerator*UNIT+eh.denominator-1)//eh.denominator
    endpoint=128*delta*3**b*((delta/2)**e*2**a+(delta/2)**a*2**e)
    Glo=Q(Gticks,UNIT)
    Ehi=Q(Eticks,UNIT)+endpoint+negative(e,a,b)
    assert Ehi<Q(limit,1000)*Glo
    print((a,e,b,p),"relative error <",str(limit)+"/1000")

certify(16,128,3,16,2048,602)
certify(17,129,3,16,2048,579)
certify(8,8,128,24,4096,750)
print("PASS")