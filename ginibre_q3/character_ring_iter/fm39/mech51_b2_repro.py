import argparse
from collections import defaultdict
from math import comb
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap=argparse.ArgumentParser(
    description="FM-MECH51: b=2 certificates and radial cutoff audit.")
ap.add_argument("--box",type=int,default=12,
                help="size of bounded Catalan bridge check")
args=ap.parse_args()
n,j,k,d,r,D=s.symbols("n j k d r D")

weights=defaultdict(int)
for ii in range(3):
    for h1,h2,sgn in ((0,0,1),(-1,1,-1),(1,1,-1),(0,2,1)):
        for h in range(ii+1):
            for ell in range(ii+1):
                pair=tuple(sorted((2+ii+h1-2*h,2+ii+h2-2*ell)))
                weights[pair]+=(
                    sgn*comb(2,ii)*2**(2-ii)*comb(ii,h)*comb(ii,ell))
weights={p:c for p,c in weights.items() if c}

den={-1:n+1,0:s.Integer(1),1:s.Integer(1),2:j+2}
num={-1:d-(j+1)*r,0:s.Integer(1),1:r,2:d*r-n}
for h in range(2,6):
    den[h+1]=den[h]*(j+h+1)
    num[h+1]=s.expand(d*num[h]-(n-h+1)*(j+h)*num[h-1])
H=(n+1)*s.prod((j+h)**2 for h in range(2,6))*(j+6)
Q=s.Poly(s.expand(sum(
    c*s.cancel(H/(den[h]*den[ell]))*num[h]*num[ell]
    for (h,ell),c in weights.items())),r)
A,B,C=[s.expand(Q.nth(h).subs(j,n+k)) for h in range(3)]
common=s.prod((n+k+h)**2 for h in range(2,6))
disc=s.cancel((B*B-4*A*C)/common)
pd=s.Poly(disc,d)
assert all(h[0]%2==0 for h,c in pd.terms())
J=s.Poly(sum(c*D**(h[0]//2) for h,c in pd.terms()),D)
assert J.degree()==5
print("quadratic and reduced quintic: PASS",flush=True)

R,u,v,w=ring("u,v,w",QQ)
def ev(expr,variables,images):
    pp=s.Poly(expr,*variables)
    powers=[[a**h for h in range(pp.degree(x)+1)]
            for x,a in zip(variables,images)]
    ans=R.zero
    for mon,c in pp.terms():
        term=R(QQ.from_sympy(c))
        for a,h in zip(powers,mon):
            term*=a[h]
        ans+=term
    return ans
def check(p,label):
    cs=list(p.values())
    assert cs and all(c>0 for c in cs),label
    if not label.startswith("outer"):
        assert p.get((0,0,0),QQ.zero)>0,label
    print(label,len(cs),min(cs),flush=True)
def abc(nn,kk,dd):
    return [ev(p,(n,k,d),(nn,kk,dd)) for p in (A,B,C)]

nn=u+3
kk=v+3
NN=2*nn+kk
pa=s.Poly(A,d)
ar=[ev(pa.nth(2*h),(n,k),(nn,kk))*NN**(2*h) for h in range(4)]
for h in range(4):
    check(sum((QQ(comb(h,i),comb(3,i))*ar[i]
               for i in range(h+1)),R.zero),
          "A Bernstein "+str(h))

jc=[ev(J.nth(h),(n,k),(nn,kk)) for h in range(6)]
O=(nn-1)*(nn+kk+2)
n0=2*(kk+1)**2
for h in range(6):
    p=-sum((QQ(comb(h,i),comb(5,i))*jc[i]*(4*O)**i
            for i in range(h+1)),R.zero)
    if h==5:
        p+=5184*nn**9*(nn-n0)
    check(p,"coverage low "+str(h))

nn=2*(kk+1)**2+u
NN=2*nn+kk
L=3*(kk+1)**2
jc=[ev(J.nth(h),(n,k),(nn,kk)) for h in range(6)]
dc=NN**2*(L-2)**2
for h in range(6):
    p=-sum((QQ(comb(h,i),comb(5,i))*jc[i]*dc**i*L**(10-2*i)
            for i in range(h+1)),R.zero)
    check(p,"coverage high "+str(h))

# Small folded levels t=0,1.
nn=u+3
kk=v+3
NN=2*nn+kk
jj=nn+kk
for tt in (0,1):
    aa,bb,cc=abc(nn,kk,NN-2*tt)
    K,KP=(R.one,R.one) if tt==0 else (kk,kk+2)
    X=(jj+1)*K
    Y=nn*KP
    check(aa*X**2+bb*X*Y+cc*Y**2,"small t="+str(tt))

# Central intervals.
for name,tt,odd in (
        ("even",u+4,False),("odd",u+3,True),("t2",R(2),False)):
    kk=v+3
    NN=3*tt*(kk+1)**2+w
    nn=(NN-kk)/2
    aa,bb,cc=abc(nn,kk,NN-2*tt)
    if odd:
        den=3*NN*kk
        lo=(kk+2)*(3*NN-(4*tt+2)*(kk+1))
        hi=3*NN*(kk+2)
        values=(
            aa*den**2+bb*lo*den+cc*lo**2,
            aa*den**2+bb*(lo+hi)*den/2+cc*lo*hi,
            aa*den**2+bb*hi*den+cc*hi**2)
    else:
        values=(aa,aa+bb/2,aa+bb+cc)
    for h,p in enumerate(values):
        check(p,"central "+name+" "+str(h))

# Outer parametrization and degree-three Bernstein coefficients.
tt=u+2
T=2*tt+1
L=4*T
x=v
y=v+w
nnum=L+x*y
knum=(2*tt-2)*L+2*T*(x+y)
dnum=2*x*y+2*T*(x+y)
powers=[[p**h for h in range(11)]
        for p in (nnum,knum,dnum,L)]
def lift(poly):
    pp=s.Poly(poly,n,k,d)
    assert pp.total_degree()<=10
    ans=R.zero
    for mon,c in pp.terms():
        term=R(QQ.from_sympy(c))
        for a,h in zip(powers,mon):
            term*=a[h]
        ans+=term*powers[3][10-sum(mon)]
    return ans
aa,bb,cc=map(lift,(A,B,C))
rn=nnum
rd=y*(2*T+x)
for h in range(4):
    p=(aa*rd**2+QQ(h,3)*bb*rn*rd
       +QQ(h*(h-1),6)*cc*rn**2)
    check(p,"outer degree-3 "+str(h))

# Independent Catalan comparison.
def moment(h):
    return 0 if h<0 or h%2 else comb(h,h//2)//(h//2+1)
def row(a,e):
    out=[1]
    for sign,count in ((1,a),(-1,e)):
        for _ in range(count):
            out=[(out[j] if j<len(out) else 0)
                 +sign*(out[j-1] if j else 0)
                 for j in range(len(out)+1)]
    return out
def half(cc,p,b):
    N=len(cc)-1
    if (N+p)%2:
        return 0
    total=0
    M=(N+p)//2
    for ii in range(b+1):
        get=lambda h:sum(
            comb(ii,j)*(cc[h-2*j] if 0<=h-2*j<=N else 0)
            for j in range(ii+1))
        det=lambda h:get(h)**2-get(h-1)*get(h+1)
        total+=comb(b,ii)*2**(b-ii)*(det(M+ii)-det(M+ii+1))
    return total
def direct(cc,p,b):
    N=len(cc)-1
    out=0
    for j,c in enumerate(cc):
        for ell in range(p//2+1):
            up=(-1)**ell*comb(p-ell,ell)
            for h in range(b+1):
                for h2 in range(b-h+1):
                    const=comb(b,h)*comb(b-h,h2)*(-2)**(b-h-h2)
                    out+=(c*up*const
                          *moment(N-j+p-2*ell+2*h)*moment(j+2*h2))
    return out

count=0
for e in range(1,args.box+1):
    for a in range(args.box+1):
        cc=row(a,e)
        N=a+e
        for p in range(7,args.box+1):
            if (N+p)%2 or N+4-p<6:
                continue
            jj=(N+p-4)//2
            val=half(cc,p,2)
            get=lambda h:cc[h] if 0<=h<=N else 0
            assert val==sum(
                c*get(jj+h)*get(jj+ell)
                for (h,ell),c in weights.items())
            assert val==direct(cc,p,2)>=0
            assert val==half([(-1)**h*c for h,c in enumerate(cc)],p,2)
            count+=1
print("b=2 Catalan bridges",count,flush=True)

# Exact radial error checks at q=1.
z=s.symbols("z")
count=0
for NN in range(2,7):
    for b in range(2,9):
        p=s.Poly(s.expand(z**NN*(1-z)*(8*z-2)**b),z)
        integral=lambda pp,l,h:sum(
            c*(h**(j+1)-l**(j+1))/(j+1)
            for (j,),c in pp.terms())
        Jval=integral(p,s.Integer(0),s.Integer(1))
        err=s.Poly(p.as_expr()*(1-z),z)
        low=integral(err,s.Integer(0),s.Rational(1,4))
        high=integral(err,s.Rational(1,4),s.Integer(1))
        abs_error=high+(-1)**b*low
        assert Jval>0 and (NN+b-2)*abs_error<=4*Jval
        count+=1
print("radial total-variation checks",count,flush=True)

# Symbolic audit of integration by parts and reflection.
z,q,nu,bb=s.symbols("z q nu bb")
beta=4*(1+q*q)
f=beta*z-2
logder=nu/z-s.Rational(1,2)/(1-z)-q*q/(2*(1-q*q*z))
lhs=-1+(1-z)*(logder+bb*beta/f)
rhs=(nu*(1-z)/z+bb*beta*(1-z)/f-s.Rational(3,2)
     -q*q*(1-z)/(2*(1-q*q*z)))
assert s.cancel(lhs-rhs)==0
zm,zp,c=s.symbols("zm zp c")
assert s.expand(
    zp*(1-c*zp)-zm*(1-c*zm)
    -(zp-zm)*(1-c*(zp+zm)))==0
print("radial IBP and reflection identities: PASS",flush=True)
print("PASS",flush=True)
