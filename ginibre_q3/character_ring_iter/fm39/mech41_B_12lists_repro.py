"""FM-MECH41 (astra_max_ceres): B for every {1,2}-list via a positive Walsh realization of the semicircle (quantile series);
consumer checks phi_r(h_1^a hat S_2^k); the canonical extension fails at (1,1,1,3)."""
import argparse
import ast
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, lcm
from pathlib import Path

argparse.ArgumentParser(
    description="FM-MECH41: quantile coefficients, charge factors, callers, obstruction"
).parse_args()

# Rational coefficients of the inverse of integral_0^v sqrt(1-t^2/4) dt.
def pmul(a,b,N):
    out=[F(0)]*N
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<N:
                out[i+j]+=x*y
    return out

def ppow(a,k,N):
    out=[F(1)]+[F(0)]*(N-1)
    for _ in range(k):
        out=pmul(out,a,N)
    return out

K=10
rho=[F(0)]+[
    F(comb(2*j,j),16**j*(2*j-1)*(2*j+1)) for j in range(1,K+1)]
coeff=[]
for k in range(K+1):
    coeff.append(sum(F(comb(2*k+ell,ell))*ppow(rho,ell,k+1)[k]
                     for ell in range(k+1))/F(2*k+1))
assert min(coeff)>0
N=2*K+2
h=[F(0)]*N
for k,v in enumerate(coeff):
    h[2*k+1]=v
back=h[:]
for j in range(1,K+1):
    term=ppow(h,2*j+1,N)
    back=[x-rho[j]*y for x,y in zip(back,term)]
assert back==[F(0),F(1)]+[F(0)]*(N-2)
print("Inverse coefficients:",*[str(x) for x in coeff[:6]])
print("Formal inverse checked through degree",N-1)

# A rational finite Boolean model checks the charge-factor identity.
# It is separate from the actual semicircle caller checks below.
G=8
def conv(a,b):
    out=[0]*G
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i^j]+=x*y
    return out

U=[F(0)]*G
U[0],U[2],U[4]=F(1,2),F(1,4),F(1,8)
U3=conv(conv(U,U),U)
v=[x+coeff[1]*y for x,y in zip(U,U3)]
f1=[v[i^1] for i in range(G)]
f2=conv(v,v)
f2[0]=F(0)
symbols={}
for n,p in ((1,f1),(2,f2)):
    den=lcm(*(x.denominator for x in p))
    symbols[n]=[int(x*den) for x in p]
    assert min(symbols[n])>=0 and symbols[n][0]==0

def rowspace(rows):
    piv={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in piv:
                x^=piv[p]
            else:
                piv[p]=x
                break
    for p in sorted(piv):
        for q in piv:
            if q>p and (piv[q]>>p)&1:
                piv[q]^=piv[p]
    return tuple(piv[p] for p in sorted(piv))

@lru_cache(None)
def finite_moment(labels):
    z=[1]+[0]*(G-1)
    for n in labels:
        z=conv(z,symbols[n])
    return z[0]

words=assignments=entries=0
for L in range(1,8):
    for n1 in range(L+1):
        labels=(1,)*n1+(2,)*(L-n1)
        supports=[[g for g,x in enumerate(symbols[n]) if x] for n in labels]
        weights={}
        for prefix in product(*supports[:-1]):
            last=0
            weight=1
            for n,g in zip(labels,prefix):
                last^=g
                weight*=symbols[n][g]
            weight*=symbols[labels[-1]][last]
            if not weight:
                continue
            charges=prefix+(last,)
            rows=rowspace([
                sum(((g>>b)&1)<<i for i,g in enumerate(charges))
                for b in range(3)])
            weights[rows]=weights.get(rows,0)+weight
            assignments+=1
        for S in range(1<<L):
            lhs=sum(w for rows,w in weights.items()
                    if all((S&r).bit_count()%2==0 for r in rows))
            left=tuple(n for i,n in enumerate(labels) if (S>>i)&1)
            right=tuple(n for i,n in enumerate(labels) if not ((S>>i)&1))
            assert lhs==finite_moment(left)*finite_moment(right)
            entries+=1
        words+=1
assert (words,assignments,entries)==(35,8931,1792)
print("Finite charge identity:",words,assignments,entries)

@lru_cache(None)
def cat(n):
    return comb(2*n,n)//(n+1)

@lru_cache(None)
def cm(p,k):
    if p%2:
        return 0
    return sum((-1)**(k-j)*comb(k,j)*cat(p//2+j) for j in range(k+1))

@lru_cache(None)
def fusion_moment(p,k):
    row={0:1}
    for n in (1,)*p+(2,)*k:
        nxt={}
        for j,v in row.items():
            for ell in range(abs(j-n),j+n+1,2):
                nxt[ell]=nxt.get(ell,0)+v
        row=nxt
    return row.get(0,0)

checks=0
for p in range(19):
    for k in range(11):
        assert cm(p,k)==fusion_moment(p,k)
        checks+=1
assert checks==209
print("Catalan/fusion moment checks:",checks)

def cb(n,k):
    return comb(n,k) if 0<=k<=n else 0

@lru_cache(None)
def kraw(n,t,s):
    return sum((-1)**j*cb(t,j)*cb(n-t,s-j) for j in range(s+1))

def fourier(p,k,t1,t2):
    return sum(kraw(p,t1,s)*kraw(k,t2,j)*cm(s,j)*cm(p-s,k-j)
               for s in range(p+1) for j in range(k+1))

checks=0
for L in range(1,13):
    for p in range(L+1):
        k=L-p
        for t1 in range(p+1):
            for t2 in range(k+1):
                z=fourier(p,k,t1,t2)
                assert z>=0 and ((t1+t2)%2==0 or z==0)
                checks+=1
assert checks==1819
print("Actual mixed-label Fourier checks:",checks)

def semicircle(n):
    return 0 if n%2 else cat(n//2)

def direct(r,a,k):
    ans=0
    for i in range(2*r+1):
        for j in range(a+1):
            for b in range(k+1):
                for c in range(k-b+1):
                    weight=((-1)**i*comb(2*r,i)*comb(a,j)
                            *comb(k,b)*comb(k-b,c)*(-2)**(k-b-c))
                    ans+=weight*semicircle(2*r-i+a-j+2*b)*semicircle(i+j+2*c)
    return ans

checks=0
for r in range(1,6):
    for a in range(7):
        for k in range(7):
            z=fourier(2*r+a,k,2*r,0)
            assert z==direct(r,a,k) and z>=0
            checks+=1
assert checks==245
print("General-sector definition checks:",checks)
checks=0
for r in range(1,9):
    for k in range(9):
        z=fourier(0,2*r+k,0,2*r)
        assert z==direct(r,2*r,k) and z>=0
        checks+=1
assert checks==72
print("Repeated-2 consumer checks:",checks)
print("Sample phi values:",
      [F(direct(r,a,k),2) for r,a,k in ((1,2,3),(2,4,5),(4,8,7))])

# Exact negative bound for the grouped E_4 coefficient.
# Integral x^4(1-x)^4/(1+x^2) = 22/7-pi > 0.
quot=[F(4),F(0),F(-4),F(0),F(5),F(-4),F(1)]
numerator=pmul(quot,[F(1),F(0),F(1)],9)
numerator[0]-=4
assert numerator==[F(0)]*4+[F(1),F(-4),F(6),F(-4),F(1)]
assert sum(v/F(i+1) for i,v in enumerate(quot))==F(22,7)
assert F(22,7)**2<10
I1=F(8,3)
I3=I1*F(8,5)
assert I3-2*I1==F(-16,15)
assert cat(3)-4*cat(2)+4*cat(1)==1
xlo=F(32,45)
A=F(2,5)*xlo*xlo
margin=(A-F(1,25))**2-(1-xlo)**3
assert A>F(1,25) and margin==F(227824,102515625)>0
print("Grouped E4 coefficient < -1/25; squared margin:",margin)

# Check the actual failed-rule table, including every q coefficient.
path=Path("ginibre_q3/character_ring_iter/fm39/sec113_hACq_L8_repro.py")
names={"trim","padd","psub","pmul","pshift","qbinom",
       "qfactorial","linearization","moment"}
tree=ast.parse(path.read_text())
env={"lru_cache":lru_cache}
exec(compile(ast.Module(body=[
    node for node in tree.body
    if isinstance(node,ast.FunctionDef) and node.name in names
],type_ignores=[]),str(path),"exec"),env)
labels=(1,1,1,3)
for S in range(16):
    left=tuple(n for i,n in enumerate(labels) if (S>>i)&1)
    right=tuple(n for i,n in enumerate(labels) if not ((S>>i)&1))
    value=env["pmul"](env["moment"](left),env["moment"](right))
    assert value==((1,2,2,1) if S in (0,15) else (0,))
print("Failed-rule table: f_q = (1+2q+2q^2+q^3) 1_<1111>")

# The three-plane relation repairs this smallest failure.
planes=[{0,15,p,15^p} for p in (3,5,9)]
for S in range(16):
    assert sum(S in H for H in planes)==int(S.bit_count()%2==0)+2*int(S in (0,15))
print("Three-plane repair: 16 exact entries")
print("All assertions passed.")