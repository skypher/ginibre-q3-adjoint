"""FM-MECH70: exact certificates and bounded verification. No file writes."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, lcm, prod
import sympy as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--part", choices=("all", "symbolic", "vet"), default="all")
args = ap.parse_args()
# e: (anchor cap, checkpoint, kappa lower bound, exp lower bound, r, w)
CFG = {
    4: (F(8,15), F(7,8), F(19,2), 96, F(1,78), F(2,9)),
    5: (F(3,5), F(7,8), F(23,2), 106, F(1,90), F(2,9)),
    6: (F(13,20), F(9,10), F(27,2), 185, F(1,160), F(9,50)),
}
EXPECTED = {4:(29,(6,9),10), 5:(13,(3,4),4), 6:(37,(10,9),10)}

def bernstein(poly):
    ds = tuple(poly.degree(z) for z in poly.gens)
    a = {idx:F(co) for idx,co in poly.terms()}
    for ax,d in enumerate(ds):
        b = {}
        for idx,co in a.items():
            h = idx[ax]
            for i in range(h,d+1):
                ix = list(idx); ix[ax] = i; ix = tuple(ix)
                b[ix] = b.get(ix,F(0)) + co*F(comb(i,h),comb(d,h))
        a = b
    den = lcm(*(co.denominator for co in a.values()))
    keys = product(*(range(d+1) for d in ds))
    return ds,[int(a.get(ix,F(0))*den) for ix in keys]

def split(ds,a,ax):
    d = ds[ax]; stride = prod(h+1 for h in ds[ax+1:])
    block = (d+1)*stride
    left = [0]*len(a); right = [0]*len(a)
    for start in range(0,len(a),block):
        for off in range(stride):
            row = [a[start+off+i*stride] for i in range(d+1)]
            for i in range(d+1):
                left[start+off+i*stride] = (1<<(d-i))*sum(
                    comb(i,h)*row[h] for h in range(i+1))
                right[start+off+i*stride] = (1<<i)*sum(
                    comb(d-i,h-i)*row[h] for h in range(i,d+1))
    return left,right

def symbolic():
    sig,V,W,X,A = S.symbols("sig V W X A")
    K = [S.Integer(1),X]
    for l in range(1,6):
        K.append(S.expand(X*K[-1]-l*(sig-l+1)*K[-2]))
    ex = lambda z: sum(z**h/factorial(h) for h in range(41))
    for e,(xm,gm,k,cut,rb,wb) in CFG.items():
        f = S.Poly(K[e]/(X if e%2 else 1),X)
        ff = sum(co*A**(power[0]//2) for power,co in f.terms())
        cap = xm*xm*4*(e+1)*(sig-e-1)
        root = S.Poly(S.expand(ff.subs(A,cap+W).subs(sig,256+V)),V,W)
        assert all(co>0 for co in root.coeffs())
        assert F(2*(e+1)*(255-e),257)>=k
        zm = gm+F(1,32)
        assert all(k*z*(1-z*z)>1 for z in (gm,zm))
        assert ex(k*xm*xm)>1+xm
        assert ex(k*(gm*gm-xm*xm))>cut
        rp = (1+gm)/(1+xm)
        wp = rp/((1-xm)*(1-gm))
        assert rp/cut<rb and wp/cut<wb
        margin = (1-rb)**2-4*F(33,32)**2*wb
        assert margin>0 and wb*(1+1/k)<1

        u,v,w = S.symbols("u v w")
        x = S.Rational(xm.numerator,xm.denominator)*u
        y = x+(S.Rational(zm.numerator,zm.denominator)-x)*v
        t = w/256
        cb = 4*(e+1)*(1-(e+1)*t)
        def kp(A):
            p = [S.Integer(1),S.Integer(1)]
            for l in range(1,e):
                p.append(S.expand((A*p[-1] if l%2 else p[-1])
                                   -l*(1-(l-1)*t)*p[-2]))
            return p
        p,q = kp(cb*x*x),kp(cb*y*y)
        G = S.Poly(S.expand(sum(
            S.Rational(factorial(e),factorial(l))*
            prod(1-h*t for h in range(l,e))*p[l]*q[l]
            for l in range(e%2,e+1,2))),u,v,w)
        kk = S.Rational(k.numerator,k.denominator)
        z = kk*(y*y-x*x)
        P = S.Poly(S.expand(
            (1-y)*(1-x*x)*sum(z**h/factorial(h) for h in range(9))
            -(1+y)*(1-x*x+1/kk)),u,v)
        gd,ga = bernstein(G); pd,pa = bernstein(P)
        stack = [(ga,pa,0)]
        leaves = [0,0]; depth = nodes = 0
        while stack:
            g,p,dep = stack.pop()
            nodes += 1
            if min(g)>0:
                leaves[0] += 1; depth=max(depth,dep); continue
            if min(p)>=0:
                leaves[1] += 1; depth=max(depth,dep); continue
            assert dep<20,("unresolved polynomial box",e,dep)
            ax = dep%2
            gl,gr = split(gd,g,ax); pl,pr = split(pd,p,ax)
            stack.append((gr,pr,dep+1))
            stack.append((gl,pl,dep+1))
        assert (nodes,tuple(leaves),depth)==EXPECTED[e]
        print("E",e,"ROOT TERMS",len(root.terms()),"J MARGIN",margin,
              "COVER",nodes,tuple(leaves),depth,flush=True)

def row(a,e):
    N=a+e; sig=N+2; de=a-e; C=4*(a+1)*(e+1)
    c=[1]
    for k in range(N):
        h,rem=divmod(de*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
        assert rem==0
        c.append(h)
    v=lambda k:c[k] if 0<=k<=N else 0
    B=[v(k-1)+v(k+1) for k in range(N+2)]
    D=[v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H=[sig*(v(k)**2+v(k-1)**2)-2*de*v(k)*v(k-1)
       for k in range(N+2)]
    beta=[comb(sig,k+1) for k in range(N+2)]
    return N,sig,C,v,B,D,H,beta

def be(C,sig,j,i,beta):
    X=2*j-sig+2; Y=2*i-sig+2
    P=C-X*X; V=P+2*(sig-X)
    h=P*beta[j]-V*beta[i]
    z=Y*(P*beta[j]+V*beta[i])
    return h>=0 and C*h*h>=z*z

@lru_cache(None)
def moment(power,label):
    if (power+label)%2:return 0
    return sum((-1)**h*comb(label-h,h)*
               (comb(power+label-2*h,(power+label-2*h)//2)//
                ((power+label-2*h)//2+1))
               for h in range(label//2+1))

@lru_cache(None)
def kernel(e,a,p,q):
    return sum((-1)**h*comb(e,h)*comb(a,k)*
               moment(e+a-h-k,p)*moment(h+k,q)
               for h in range(e+1) for k in range(a+1))

def vet():
    rows=pairs=0; least=None; local_counts=[]
    # The remaining fixed box after MECH61's omega<64 theorem.
    for e in (4,5,6):
        local=0
        for a in range(e+2,254-e):
            N,sig,C,v,B,D,H,beta=row(a,e)
            if C<4096:continue
            rows+=1
            for j in range((N+1)//2,N-3):
                for i in range(j+4,N+2):
                    slack=D[j]-D[i]-abs(v(j)*B[i]-B[j]*v(i))
                    assert slack>=0,(a,e,j,i)
                    local+=1; pairs+=1
                    item=slack,(a,e,j,i)
                    if least is None or item<least:least=item
        local_counts.append(local)
    assert local_counts==[297252,442586,515848]
    assert (rows,pairs)==(227,1255686)
    assert least==(48132946518867,(146,6,148,152))
    print("FIXED BOX",rows,pairs,local_counts,"MINIMUM",least,flush=True)

    bridges=0
    for e in (4,5,6):
        for a in range(e+2,25):
            N,sig,C,v,B,D,H,beta=row(a,e)
            def K(X):
                p=[1,X]
                for l in range(1,e):
                    p.append(X*p[-1]-l*(sig-l+1)*p[-2])
                return p
            ks={j:K(2*j-N) for j in range((N+1)//2,N+2)}
            eta=F(factorial(e),4*prod(sig-h for h in range(e+2)))
            nu=[factorial(l)*prod(sig-h for h in range(l))
                for l in range(e+1)]
            for j in range(N//2+1,N):
                for i in range(j+1,N+2):
                    X,Y=2*j-N,2*i-N
                    G=sum(F(ks[j][l]*ks[i][l],nu[l])
                          for l in range(e%2,e+1,2))
                    W=v(j)*B[i]-B[j]*v(i)
                    assert W==eta*(Y*Y-X*X)*beta[j]*beta[i]*G
                    bridges+=1
    assert bridges==3380
    print("RF BRIDGES",bridges,flush=True)

    direct=0
    for e in (4,5,6):
        for a in range(e+2,14):
            N,sig,C,v,B,D,H,beta=row(a,e)
            for j in range(N//2+1,N-3):
                for i in range(j+4,N+2):
                    p,q=i+j-N-1,i-j-1
                    W=v(j)*B[i]-B[j]*v(i)
                    assert kernel(e,a,p,q)==W
                    assert sum(kernel(e,a,l,0)
                               for l in range(p-q,p+q+1,2))==D[j]-D[i]
                    direct+=1
    assert direct==273
    print("CATALAN BRIDGES",direct,flush=True)

    sample=cross=tail=0
    for e,(xm,gm,*_) in CFG.items():
        for a in (255,1000,4000):
            N,sig,C,v,B,D,H,beta=row(a,e)
            for j in range(N//2+1,N-3):
                X=2*j-N
                if X*X>=xm*xm*C:break
                first=False
                for i in range(j+4,N+2):
                    Y=2*i-N
                    if Y*Y>F(11,10)**2*C:break
                    W=v(j)*B[i]-B[j]*v(i)
                    drop=D[j]-D[i]
                    assert drop>=abs(W),(a,e,j,i)
                    sample+=1
                    if not first and W<0 and Y*Y<gm*gm*C:
                        assert be(C,sig,j,i,beta),(a,e,j,i)
                        first=True; cross+=1
                    if Y*Y>=gm*gm*C:
                        assert C*(sig+Y)**2*drop**2>=4*Y*Y*H[j]*H[i]
                        tail+=1
    assert (sample,cross,tail)==(43383,424,12137)
    print("LARGE SAMPLES",sample,"FIRST BE",cross,"TAIL J",tail)

if args.part in ("all","symbolic"):symbolic()
if args.part in ("all","vet"):vet()
print("PASS")

