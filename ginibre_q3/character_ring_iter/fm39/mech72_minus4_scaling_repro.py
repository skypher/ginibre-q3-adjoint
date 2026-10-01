import argparse
from fractions import Fraction as Q
from math import comb,factorial,isqrt

argparse.ArgumentParser(
    description='FM-MECH72 exact Gaussian transfer verifier'
).parse_args()

def add(a,b):
    return a[0]+b[0],a[1]+b[1]
def scale(a,c):
    return min(a[0]*c,a[1]*c),max(a[0]*c,a[1]*c)
def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)

def exneg(x):
    assert 0<=x<=1
    lo=sum(((-x)**j/Q(factorial(j))
            for j in range(40)),Q(0))
    return lo,lo+x**40/Q(factorial(40))

def atan_small(d):
    lo=sum((Q((-1)**j,(2*j+1)*d**(2*j+1))
            for j in range(24)),Q(0))
    return lo,lo+Q(1,49*d**49)

# Machin's identity and an integer square-root enclosure.
p5,p239=atan_small(5),atan_small(239)
pi=(16*p5[0]-4*p239[1],16*p5[1]-4*p239[0])
unit=10**40
sq=isqrt(3*unit*unit)
root=(Q(sq,unit),Q(sq+1,unit))
c0=(Q(125,512)/(pi[1]*root[1]),
    Q(125,512)/(pi[0]*root[0]))

def H(a,b):
    z=4*a*b
    lo=sum((z**j/Q(factorial(2*j+1))
            for j in range(20)),Q(0))
    tail=z**20/Q(factorial(41))
    ratio=z/Q(42*43)
    assert ratio<1
    return mul(exneg(a+b),(lo,lo+tail/(1-ratio)))

def va(n):
    # For odd n this is U_n(sqrt(3/2))/((n+1)*sqrt(3/2)).
    return sum((
        Q((-1)**j*comb(n-j,j))*Q(3,2)**(n//2-j)
        for j in range(n//2+1)
    ),Q(0))/(n+1)

def gaussian(n,h,m=None,eps2=1):
    r=Q(h,2)
    a=Q(5*(n+1)**2,128)/r
    if m is None:
        return mul(c0,add(exneg(a),(((-1)**h)*va(n),)*2))
    b=Q(5*(m+1)**2,128)/r
    eps1=(-1)**h*eps2
    g=H(a,b)
    if n%2:
        prod=eps1*eps2*Q(3,2)*va(n)*va(m)
        g=add(g,(prod,prod))
    else:
        g=add(g,scale(exneg(a),eps2*va(m)))
        g=add(g,scale(exneg(b),eps1*va(n)))
        prod=eps1*eps2*va(n)*va(m)
        g=add(g,(prod,prod))
    return mul(c0,g)

def exact_values(n,h,m):
    # Character coefficients of U_4^j, indexed by half the label.
    row=[1]
    z=[];a=[];b=[];ab=[]
    for j in range(h+1):
        z.append(row[0])
        a.append(row[n//2]
                 if n%2==0 and n//2<len(row) else 0)
        b.append(row[m//2]
                 if m%2==0 and m//2<len(row) else 0)
        ab.append(sum(row[abs(n-m)//2:(n+m)//2+1]))
        if j==h:
            break
        nxt=[0]*(len(row)+2)
        for i,v in enumerate(row):
            for ell in range(abs(i-2),i+3):
                nxt[ell]+=v
        row=nxt
    one=2*sum(
        (-1)**j*comb(h,j)*a[h-j]*z[j]
        for j in range(h+1)
    )
    common=2*sum(
        (-1)**j*comb(h,j)*ab[h-j]*z[j]
        for j in range(h+1)
    )
    cross=2*sum(
        (-1)**j*comb(h,j)*a[h-j]*b[j]
        for j in range(h+1)
    )
    return one,common+cross,common-cross

# Constants in the local Taylor estimates and global tails.
assert Q(7,400)<Q(1,50)
assert Q(50,49)*Q(32,625)*441<24
assert Q(50,49)*Q(32,625)*49<3
assert Q(8,25)*23+24<32
assert Q(8,25)*6+Q(3,16)<3
assert 23*512>6400
R=2**36
err=(Q(1000,R)+Q(24,2**18)
     +Q(8*512**2,R**2)+Q(480*800**5,R**3))
assert err<Q(1,10000)
assert Q(125,4096)*Q(5,56)**2>Q(1,8192)
assert Q(125,4096)*Q(19,72)>Q(1,8192)
assert Q(8192,10000)==Q(512,625)<1
print('Analytic constants and relative-error cutoff: PASS')

count=0
for n in (5,6,10,14):
    m=n+2
    M=n+1
    for h in (2*M*M,2*M*M+1):
        r=Q(h,2)
        one,plus,minus=exact_values(n,h,m)
        data=[]
        if n%2==0:
            data.append((one,M,gaussian(n,h)))
        data += [
            (F,M*(m+1),gaussian(n,h,m,eps))
            for F,eps in ((plus,1),(minus,-1))
        ]
        for F,dim,G in data:
            Z=r*r*Q(F*4**h,dim*25**h)
            assert F>0 and G[0]>0
            assert Q(97,100)*G[1]<Z<Q(103,100)*G[0]
            count+=1
print('Exact edge-model ratio brackets (97/100,103/100):',count)

def C(n,k):
    return comb(n,k) if 0<=k<=n else 0

def ceres_value(p):
    b=p*p
    T=p+2
    row=[C(p,j)-2*C(p,j-1)+C(p,j-2)
         for j in range(T+1)]
    cats=[comb(2*j,j)//(j+1)
          for j in range(b+T//2+1)]
    ballot={
        m:(C(m,(m-p)//2)-C(m,(m-p)//2-1)
           if m>=p and (m-p)%2==0 else 0)
        for m in range(T+2*b+1)
    }
    return sum(
        row[j]*C(b,l)*(-2)**(b-l)
        *sum(C(l,h)*cats[j//2+h]*ballot[T-j+2*(l-h)]
             for h in range(l+1))
        for j in range(0,T+1,2)
        for l in range(b+1)
    )

for p,digits in ((8,13),(10,21),(12,31)):
    F=ceres_value(p)
    b=p*p
    d=b+1
    T=p+2
    lower=Q((b-T//2)**d*comb(2*d,d),
            factorial(d)*(d+1))
    assert 0<Q(F)/lower<Q(1,10**digits)
    a=Q(p,4)+Q(2*b,3)
    t=Q((p+1)**2)/(4*a)
    ep=exneg(t)
    P=t*t-2*t+3
    G=(ep[0]*P/pi[1],ep[1]*P/pi[0])
    Z=2*a**5*Q(F,4**p*6**b*(p+1))
    assert Q(9,10)*G[1]<Z<G[0]
    print('MECH53 p',p,'distance',d,
          'Ld ratio < 10^-'+str(digits),
          'corner ratio in (9/10,1)')

# Corner obstruction: n=6, r=49, angular radius 1/7.
M=7
r=M*M
corner_ratio=Q(393216*M**3,7)*Q(8,M*M)**(2*r)
assert corner_ratio<Q(1,10**69)
assert 2*r<384*(2*M*M)
assert Q(2*M*M+50*r,M*M)<2**21
print('Corner fraction < 10^-69 at n=6,r=49: PASS')

# (1-v/4)^4(1+v), with v=y^2.
p=[Q(comb(4,j)*(-1)**j,4**j) for j in range(5)]
critical=[
    (p[j] if j<5 else 0)+(p[j-1] if j else 0)
    for j in range(6)
]
assert critical==[
    1,0,Q(-5,8),Q(5,16),Q(-15,256),Q(1,256)
]
print('Critical mixed-background quartic identity: PASS')
print('PASS')
