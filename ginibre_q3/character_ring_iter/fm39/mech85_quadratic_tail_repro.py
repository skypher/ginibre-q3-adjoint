import argparse
argparse.ArgumentParser(
    description="FM-MECH85 exact integrated Gram tests; memory only"
).parse_args()
from math import comb
from fractions import Fraction as Q
from functools import lru_cache
import sympy as S,time

def mul(a,b,cap):
 out={}
 for (i,j),v in a.items():
  for (k,l),w in b.items():
   if i+k+2*(j+l)<=cap:
    ij=(i+k,j+l);out[ij]=out.get(ij,0)+v*w
 return {m:v for m,v in out.items() if v}

def power(a,n,cap):
 ans={(0,0):1}
 while n:
  if n&1:ans=mul(ans,a,cap)
  n//=2
  if n:a=mul(a,a,cap)
 return ans

def block(n,eps,cap):
 f={(0,j):1 for j in range(min(n,cap//2)+1)}
 if n<=cap:
  for j in range(n//2+1):
   m=(n-2*j,j);f[m]=f.get(m,0)+eps*(-1)**j*comb(n-j,j)
 return {m:v for m,v in f.items() if v}

def kernel(counts,D):
 p={(0,0):1}
 for n,e,c in counts:
  p=mul(p,power(block(n,e,2*D),c,2*D),2*D)
 return S.Matrix(D+1,D+1,lambda i,j:sum(
  S.Rational(2*comb(i+j-2*h,i-h),i+j-2*h+2)*p.get((i+j-2*h,h),0)
  for h in range(min(i,j)+1)))

def pivots(A):
 A=[list(A.row(j)) for j in range(A.rows)]
 out=[]
 for j in range(len(A)):
  p=A[j][j];out.append(p)
  if p<0:return out
  if p==0:
   if any(A[j][i] for i in range(j+1,len(A))):
    return out+["zero-cross"]
   continue
  for a in range(j+1,len(A)):
   for b in range(a,len(A)):
    v=A[a][b]-A[a][j]*A[j][b]/p
    A[a][b]=A[b][a]=v
 return out

def shifted(A):
 return S.Matrix(A.rows,A.cols,lambda i,j:
  A[i,j]-(A[i-1,j-1] if min(i,j)>0 else 0))

def pointkernel(counts,D,c):
 p={(0,0):1}
 for n,e,count in counts:
  p=mul(p,power(block(n,e,2*D),count,2*D),2*D)
 return S.Matrix(D+1,D+1,lambda i,j:sum(
  comb(i+j-2*h,i-h)*c**(i+j-2*h)*p.get((i+j-2*h,h),0)
  for h in range(min(i,j)+1)))

def shifted_by(A,scale):
 return S.Matrix(A.rows,A.cols,lambda i,j:
  A[i,j]-(scale*A[i-1,j-1] if min(i,j)>0 else 0))

def positive(A):
 return all(v!="zero-cross" and v>0 for v in pivots(A))

from random import Random
c=S.symbols('c')
for n in range(1,19):
 for j in range(n+1):
  a=sum((-1)**l*comb(n-l,l)*comb(n-2*l,j-l)*c**(n-2*l)
        for l in range(min(j,n-j)+1))
  d=sum((-1)**l*comb(j,l)*comb(n-j,l)*c**(n-2*l)*(1-c*c)**l
        for l in range(min(j,n-j)+1))
  assert S.expand(a-d)==0
print("Wigner diagonal polynomial identities: 189 PASS",flush=True)

def vals(cs,D):
 p={(0,0):1}
 for n,e,count in cs:
  p=mul(p,power(block(n,e,2*D),count,2*D),2*D)
 return [sum(v*comb(i,i//2)//(i//2+1) for (i,j),v in p.items()
             if i%2==0 and i//2+j==d) for d in range(D+1)]

@lru_cache(None)
def mult(ns):
 if not ns:return 1
 if sum(ns)%2 or 2*max(ns)>sum(ns):return 0
 row={0:1}
 for n in ns:
  rr={}
  for j,a in row.items():
   for m in range(abs(j-n),j+n+1,2):
    rr[m]=rr.get(m,0)+a
  row=rr
 return row.get(0,0)

def even(word):
 L=len(word);out=0
 for mask in range(1<<L):
  A=[];B=[];sg=1
  for j,(n,e) in enumerate(word):
   if mask>>j&1:A.append(n);sg*=e
   else:B.append(n)
  out+=sg*mult(tuple(sorted(A)))*mult(tuple(sorted(B)))
 return out

rng=Random(8501)
for case in range(120):
 bg=[(rng.randrange(1,8),rng.choice((-1,1)))
     for j in range(rng.randrange(1,8))]
 D=sum(n for n,e in bg)
 pp=vals([(n,e,1) for n,e in bg],D)
 assert all(a>0 for a in pp)
 d=rng.randrange(1,13);h=rng.randrange(1,5)
 if h==1 and D<2*d:continue
 B=2*d+max(n for n,e in bg)+4
 rest=bg+[(B,1)]*h;p=sum(n for n,e in rest)-2*d
 ep=(-1)**sum(e<0 for n,e in rest)
 P=pp[:d+1]+[0]*max(0,d+1-len(pp))
 if h==1:value=P[d]
 else:value=sum(P[j]*comb(d-j+h-2,h-2) for j in range(d+1))
 assert even(rest+[(p,ep)])==2*value>0
print("120 coefficient profiles and admissible saturation bridges PASS",
      flush=True)

for r,expected in [(7,285558),(9,10658882),(11,407294418)]:
 pp=vals([(3,1,r),(3,-1,r),(4,1,1)],r)
 assert 2*sum(pp)==expected
 print("FM83 escape family",r,"value",2*sum(pp),
       "minimum P coefficient",min(pp),flush=True)

st=time.monotonic();tests=0
for d in (1,2,4,6,8,10):
 k=40*d
 profiles=[
  [(3,1,k)],
  [(3,-1,k)],
  [(3,1,k//2),(3,-1,k-k//2)],
  [(3,-1,k//3),(4,1,k//3),(7,-1,k-2*(k//3))],
  [(2,1,3),(1,1,5),(1,-1,3),(d+3,-1,k)]]
 for cs in profiles:
  for c in (S.Rational(0),S.Rational(1,3),
            S.Rational(3,5),S.Rational(1)):
   assert positive(shifted_by(pointkernel(cs,d,c),16)),(d,cs,c)
   tests+=1
  G=kernel(cs,d)
  assert positive(shifted_by(G,16)),(d,cs)
  tests+=1
 print("core shift bound d",d,"checks",tests,
       "seconds",round(time.monotonic()-st,2),flush=True)

d=S.symbols('d')
assert S.expand(4*(d+1)**3-27*d*d-(d-2)**2*(4*d+1))==0
assert 77**2*40>480**2
assert 27*432*156**2>28800*40**2*4
for d in range(1,101):
 s=d+1;N=432*s*s
 assert N*(3*N-1140*s*s)**2>=28800*(40*d)**2*s**3
print("quadratic-tail integer consumer: 100 exact cases and symbolic bound PASS",
      flush=True)

def gb(n,j):
 if j<0:return 0
 if n>=0:return comb(n,j) if j<=n else 0
 return (-1)**j*comb(j-n-1,j)

@lru_cache(None)
def moment(q,word):
 row={0:1}
 for n in (1,)*q+word:
  nxt={}
  for j,a in row.items():
   for ell in range(abs(j-n),j+n+1,2):
    nxt[ell]=nxt.get(ell,0)+a
  row=nxt
 return row.get(0,0)

def signed_ec(c,y,r):
 a,b=1,y
 if not r:return a
 for j in range(1,r):
  v=y*b-(c-j+1)*a
  assert v%(j+1)==0
  a,b=b,v//(j+1)
 return b

def deplete_value(d,cs):
 t=x=b=k=0;counts={}
 for n,e,c in cs:
  if n==1:t+=c;x+=e*c
  elif n==2:
   if e==1:b+=c
   else:t+=2*c
  else:
   k+=c
   if n<=d:
    co,si=counts.get(n,(0,0));counts[n]=(co+c,si+e*c)
 A=[1,x]
 for j in range(1,2*d):
  v=x*A[j]-(t-j+1)*A[j-1]
  assert v%(j+1)==0
  A.append(v//(j+1))
 @lru_cache(None)
 def ker(D,word):
  if D<0 or sum(word)>2*D:return 0
  w=sum(word);r=len(word);out=0
  for j in range(w%2,2*D-w+1,2):
   if not A[j]:continue
   R0=D-(w+j)//2
   for v in range(R0//2+1):
    for u in range(R0-2*v+1):
     if u+v>b:continue
     m=moment(j+2*u,word)
     R=R0-u-2*v
     rem=sum(gb(k-r+h-2,h)*gb(t-j,R-h) for h in range(R+1))
     out+=A[j]*gb(b,u)*gb(b-u,v)*m*rem
  return out
 items=sorted(counts.items())
 def rec(j,word,spent,coef):
  if j==len(items):return coef*ker(d-spent,tuple(word))
  n,(co,si)=items[j];out=0
  for a in range(min(co,(2*d-2*spent-sum(word))//n)+1):
   ea=signed_ec(co,si,a)
   if not ea:continue
   limit=(2*d-2*spent-sum(word)-n*a)//(2*(n+1))
   for beta in range(min(co-a,limit)+1):
    if sum(word)+n*a+2*(spent+(n+1)*beta)>2*d:continue
    out+=rec(j+1,word+[n]*a,spent+(n+1)*beta,
             coef*ea*(-1)**beta*comb(co-a,beta))
  return out
 return rec(0,[],0,1)

def pvalues(cs,D):
 p={(0,0):1}
 for n,e,count in cs:
  p=mul(p,power(block(n,e,2*D),count,2*D),2*D)
 return [sum(v*comb(i,i//2)//(i//2+1) for (i,j),v in p.items()
             if i%2==0 and i//2+j==d) for d in range(D+1)]

def fvalue(cs,d):
 p=pvalues(cs,d)
 return p[d]-(p[d-1] if d else 0)

def band(bg,d,h):
 J=(d-1)//3;B=pointkernel(bg,J,S.Integer(1))
 return S.Matrix(J+1,J+1,lambda i,j:sum(
  gb(h-4+r,r)*B[i-r,j-r] for r in range(min(i,j)+1)))

rng=Random(8503);st=time.monotonic()
for case in range(100):
 d=rng.randrange(1,13)
 cs=[(rng.randrange(1,d+4),rng.choice((-1,1)),rng.randrange(1,4))
     for j in range(rng.randrange(1,7))]
 assert deplete_value(d,cs)==fvalue(cs,d),(d,cs)
print("general all-order depletion: 100 exact profiles PASS",flush=True)

for d in (7,8,9,12):
 bg=[(1,1,7),(1,-1,5),(2,1,3),(3,1,2),(3,-1,1),(4,-1,2)]
 for h in (3,4,7):
  J=(d-1)//3;G=band(bg,d,h);ps=pivots(G)
  if h>=4:
   assert all(p>=gb(h+j-4,j) for j,p in enumerate(ps))
  for i in range(J+1):
   for j in range(J+1):
    c=0
    for a in (-1,1):
     for b in (-1,1):
      cs=bg+[(d-i,a,1),(d-j,b,1),(d+1,1,h-2)]
      c+=a*b*fvalue(cs,d)
    assert c==4*G[i,j],(d,h,i,j,c,G[i,j])
 print("general band Gram d",d,"PASS",flush=True)

bg=[(1,1,7),(1,-1,5),(2,1,3),(3,1,2),(3,-1,1),(4,-1,2)]
assert band(bg,7,4)==S.Matrix([[1,2,-1],[2,16,23],[-1,23,120]])
d=8;t=432*(d+1)**2
G=band([(1,1,t//2),(1,-1,t//2)],d,3)
assert G.rank()==1
print("large-N band-Gram strictness witness:",G.tolist(),flush=True)

for d in (8,10,16):
 k=40*d;cs=[(3,1,k),(1,-1,1)]
 P=pvalues(cs,d);F=P[d]-P[d-1]
 assert 16*F>=15*P[d]>0
 print("core theorem sample",d,"N",k+1,"F",F,flush=True)

d=8;s=d+1;N=432*s*s;k=40*d-1;b=3;t=N-k-2*b
cs=[(1,1,(t+1)//2),(1,-1,t//2),(2,1,b),(3,-1,k)]
F=fvalue(cs,d)
assert F>0 and N*(3*N-1140*s*s)**2>=28800*k*k*s**3
print("quadratic-tail below core-count cutoff sample F",F,flush=True)

# Shifted Krawtchouk coordinates.
d=7;J=2
bg=[(1,1,7),(1,-1,5),(2,1,3),(3,1,2),(3,-1,1),(4,-1,2)]
row=[1,2];M=12-2*3
row.append(Q(2*row[1]-M*row[0],2))
TP=S.Matrix(J+1,J+1,lambda i,j:row[i-j] if i>=j else 0)
G=band(bg,d,4)
assert TP.det()==1 and positive(TP.inv()*G*TP.inv().T)
assert moment(0,(3,3,3,3))==4
assert Q(moment(0,(3,3,3,3)),24)==Q(1,6)

def background_box(X):
 q,r=divmod(X,2)
 return (q+1)*(q+2)*(4*q+3+6*r)//6

for X in range(41):
 assert background_box(X)==sum(
  comb(X-2*b+2,2) for b in range(X//2+1))

def box(d):
 Qmax=432*(d+1)**2-1
 return sum(comb(k+2*d-5,2*d-5)*background_box(Qmax-k)
            for k in range(40*d))

assert box(8)==10274329643463376000849100761835148
print("residual box d=8:",box(8),flush=True)

def band_value(cs,d):
 J=(d-1)//3;m0=d-J
 bg=[(n,e,c) for n,e,c in cs if n<m0]
 h=sum(c for n,e,c in cs if n>=m0)
 counts=[sum(c for n,e,c in cs if n==d-i) for i in range(J+1)]
 signs=[sum(e*c for n,e,c in cs if n==d-i) for i in range(J+1)]
 pp=pvalues(bg,d)
 def base(r):
  return sum(pp[j]*gb(h+r-j-2,r-j) for j in range(r+1))
 out=Q(base(d))
 poly={(0,0):1}
 for n,e,c in bg:
  poly=mul(poly,power(block(n,e,2*d),c,2*d),2*d)
 for i,n in enumerate(range(d,d-J-1,-1)):
  term={(n-2*r,r):(-1)**r*comb(n-r,r) for r in range(n//2+1)}
  prod=mul(poly,term,2*d)
  L=0
  for (r,j),v in prod.items():
   if r%2==0:
    degree=r//2+j
    L+=v*comb(r,r//2)//(r//2+1)*gb(h+d-degree-3,d-degree)
  out+=signs[i]*L-counts[i]*base(d-n-1)
 G=band(bg,d,h)
 out+=sum(Q(signs[i]*signs[j])*G[i,j]/2
          for i in range(J+1) for j in range(J+1))
 out-=sum(Q(counts[i])*G[i,i]/2 for i in range(J+1))
 return out

rng=Random(8510)
for case in range(60):
 d=rng.randrange(3,14)
 cs=[(rng.randrange(1,d+4),rng.choice((-1,1)),rng.randrange(1,4))
     for j in range(rng.randrange(1,8))]
 assert band_value(cs,d)==fvalue(cs,d),(d,cs)
print("complete high-band quadratic identity: 60 exact profiles PASS",
      flush=True)
print("ALL CHECKS PASS",flush=True)