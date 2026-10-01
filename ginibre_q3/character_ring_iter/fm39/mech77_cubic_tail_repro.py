import argparse
from fractions import Fraction as Q
from math import comb, factorial
from functools import lru_cache
from random import Random

argparse.ArgumentParser(
    description="FM-MECH77 exact certificates; memory only"
).parse_args()

def add(P,R):
 S=P.copy()
 for m,c in R.items(): S[m]=S.get(m,Q(0))+c
 return {m:c for m,c in S.items() if c}

def scale(P,c):
 return {m:c*v for m,v in P.items() if c*v}

def timesN(P):
 R={}
 for (r,s,u),v in P.items():
  for mon in ((r+1,s,u),(r,s+1,u),(r,s,u+1)):
   R[mon]=R.get(mon,Q(0))+v
 return R

@lru_cache(None)
def alpha(d,j):
 if j<0 or j>d:return {}
 return {
  (r,s,d-j-r-s):
   Q(factorial(j+s),
     factorial(s)*factorial(d-j-r-s)*factorial(j+s+r+1))
  for r in range(d-j+1) for s in range(d-j-r+1)
 }

checks=0
for d in range(1,25):
 for j in range(d+1):
  Z=add(scale(alpha(d,j),2*(d+1)),
        scale(timesN(alpha(d-1,j)),-1))
  assert all(c>=0 for c in Z.values());checks+=1
  if j:
   Z=add(scale(alpha(d,j-1),2*(d+1)),
         scale(timesN(alpha(d,j)),-1))
   assert all(c>=0 for c in Z.values());checks+=1
print("coefficientwise norm bounds:",checks)

r=Q(1,12)
assert Q(81,4)/(1-9*r*r)<=22
assert Q(25,2)/(1-5*r*r)<=13
assert 9/(1-9*r*r)<=10
assert (29+12*r)/(1-2*r-r*r)<37
assert 1+12*r+39*r*r<Q(7,3)
assert Q(49,9)/(2*(1-Q(7,3)*r*r))<3
assert Q(99,70)**2>2
assert Q(57,4096)+Q(3,8)*Q(99,70)<Q(3,5)

for s in range(1,41):
 C=3200*s**3+760*s**2
 assert (C-380*s*s)**2-3200*C*s**3==(380*s*s)**2
 assert C>=380*s*s and C<4096*s**3

p=Q(3,5)
upper=sum(p**j/factorial(j) for j in range(9))
upper+=p**9/factorial(9)/(1-p/10)
assert upper<Q(11,6)
assert 2-(1+Q(1,144))*Q(11,6)>Q(1,8)
print("all rational constants: PASS")

def uniadd(a,b):
 c=[Q(0)]*max(len(a),len(b))
 for i,v in enumerate(a):c[i]+=v
 for i,v in enumerate(b):c[i]+=v
 return c

def umul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):c[i+j]+=v*w
 return c

def hermite_polys(d,M):
 H=[[Q(1)],[Q(0),Q(1)]]
 for j in range(1,d):
  H.append([v/(j+1) for v in
            uniadd([Q(0)]+H[j],[-M*v for v in H[j-1]])])
 return H[:d+1]

def av(d,j,t,b,k):
 return sum(c*t**r*(2*b)**s*k**u
            for (r,s,u),c in alpha(d,j).items())

def cat(j):
 return comb(2*j,j)//(j+1)

def leading(d,t,b,k,x):
 M=t-2*b;H=[Q(1),Q(x)]
 for j in range(1,2*d):
  H.append((x*H[j]-M*H[j-1])/(j+1))
 direct=sum(Q(cat(j),factorial(d-j))*(t+k)**(d-j)*H[2*j]
            for j in range(d+1))
 square=sum(av(d,j,t,b,k)*H[j]**2 for j in range(d+1))
 assert direct==square
 return square

checks=0
for d in range(13):
 for t,b,k in ((1,0,0),(0,1,0),(0,0,1),
               (7,3,2),(2,13,5),(0,4,3)):
  H=hermite_polys(2*d,t-2*b);left=[Q(0)]
  for j in range(d+1):
   left=uniadd(left,
              [av(d,j,t,b,k)*v for v in umul(H[j],H[j])])
  right=[Q(0)]
  for j in range(d+1):
   right=uniadd(right,
    [Q(cat(j),factorial(d-j))*(t+k)**(d-j)*v for v in H[2*j]])
  assert left==right;checks+=1
print("Hermite identity, full x polynomials:",checks)

def gb(n,j):
 if j<0:return 0
 if n>=0:return comb(n,j) if j<=n else 0
 return (-1)**j*comb(j-n-1,j)

def uchar(n):
 return {n-2*j:(-1)**j*comb(n-j,j) for j in range(n//2+1)}

def pdmul(a,b,cap):
 c={}
 for (i,j),v in a.items():
  for (k,l),w in b.items():
   if i+k<=cap:
    m=(i+k,j+l);c[m]=c.get(m,0)+v*w
 return {m:v for m,v in c.items() if v}

def yproduct(a,b):
 c={}
 for i,v in a.items():
  for j,w in b.items():c[i+j]=c.get(i+j,0)+v*w
 return {j:v for j,v in c.items() if v}

def block(n,eps,count,cap):
 if not count:return {(0,0):1}
 out={};yp={0:1};un=uchar(n) if n<=cap else {}
 for j in range(min(count,cap//n)+1):
  count0=count-j
  for q in range((cap-n*j)//2+1):
   u=sum((-1)**h*gb(count0,h)*
         gb(count0+q-(n+1)*h-1,q-(n+1)*h)
         for h in range(q//(n+1)+1))
   if u:
    for rr,v in yp.items():
     m=(n*j+2*q,rr)
     out[m]=out.get(m,0)+comb(count,j)*eps**j*u*v
  yp=yproduct(yp,un)
 return {m:v for m,v in out.items() if v}

def distance_value(d,counts):
 p={(0,0):1}
 for n,eps,c in counts:
  p=pdmul(p,block(n,eps,c,2*d),2*d)
 return sum(v*cat(j//2)*((i==2*d)-(i==2*d-2))
            for (i,j),v in p.items() if j%2==0)

def resources(counts):
 t=x=b=k=s3=0
 for n,eps,c in counts:
  if n==1:t+=c;x+=eps*c
  elif n==2 and eps==-1:t+=2*c
  elif n==2:b+=c
  else:k+=c;s3+=eps*c if n==3 else 0
 return t,b,k,x,s3

def criterion(d,counts):
 t,b,k,x,s3=resources(counts);N=t+2*b+k;s=d+1
 return (N>=380*s*s and
         N*(3*N-1140*s*s)**2>=28800*s3*s3*s**3)

checks=0
for d in (3,6,7,8,10,12):
 C=3200*(d+1)**3+760*(d+1)**2
 rows=[
  [(1,1,C//2),(1,-1,C-C//2),(2,1,3),(3,-1,5),(4,1,2)],
  [(1,1,7),(1,-1,5),(2,1,C//2),(3,1,3),(5,-1,2)],
  [(1,1,7),(1,-1,5),(2,1,3),(3,-1,C),(5,1,2)],
  [(1,1,5),(1,-1,3),(3,1,C//2),(3,-1,C-C//2),(d+1,1,2)]
 ]
 for cs in rows:
  t,b,k,x,s3=resources(cs)
  v=distance_value(d,cs);l=leading(d,t,b,k,x)
  assert criterion(d,cs) and 8*v>=l>0
  checks+=1
print("cubic-tail exact profiles:",checks)

checks=0
for d in (7,8,10,12,16,20):
 C=380*(d+1)**2
 rows=[
  [(1,1,C//2),(1,-1,C//2),(2,-1,3),(3,1,4),(3,-1,4)],
  [(1,1,4),(1,-1,2),(2,1,C//2),(4,-1,11),(9,1,3)],
  [(1,1,5),(1,-1,3),(3,1,C//2),(3,-1,C//2),(d+3,-1,2)]
 ]
 for cs in rows:
  t,b,k,x,s3=resources(cs)
  assert s3==0 and criterion(d,cs)
  v=distance_value(d,cs);l=leading(d,t,b,k,x)
  assert 8*v>=l>0;checks+=1
print("quadratic-tail exact profiles:",checks)

def inv_mult(ns):
 row={0:1}
 for n in ns:
  out={}
  for j,v in row.items():
   for ell in range(abs(j-n),j+n+1,2):
    out[ell]=out.get(ell,0)+v
  row=out
 return row.get(0,0)

@lru_cache(None)
def multiplicity(ns):
 if not ns:return 1
 if sum(ns)%2 or 2*max(ns)>sum(ns):return 0
 return inv_mult(ns)

def even_word(word):
 z=0;L=len(word)
 for mask in range(1<<L):
  a=[];b=[];eps=1
  for i,(n,e) in enumerate(word):
   if mask>>i&1:a.append(n);eps*=e
   else:b.append(n)
  z+=eps*multiplicity(tuple(sorted(a)))*multiplicity(tuple(sorted(b)))
 return z

rng=Random(7701);checks=0
for target,d in zip((20,40,60,80),(7,8,9,12)):
 while checks<target:
  rest=[(rng.randrange(1,13),rng.choice((-1,1)))
        for _ in range(rng.randrange(3,8))]
  n=sum(j for j,e in rest)-2*d
  if n<max(j for j,e in rest):continue
  eps=(-1)**sum(e<0 for j,e in rest)
  v=distance_value(d,[(j,e,1) for j,e in rest])
  assert 2*v==even_word(rest+[(n,eps)])
  clipped=[(min(j,d+1),e if j<=d else 1) for j,e in rest]
  nc=sum(j for j,e in clipped)-2*d
  ec=(-1)**sum(e<0 for j,e in clipped)
  assert nc>0 and 2*v==even_word(clipped+[(nc,ec)])
  split=[]
  for j,e in rest:
   split.extend([(1,1,1),(1,-1,1)]
                if (j,e)==(2,-1) else [(j,e,1)])
  assert distance_value(d,split)==v
  checks+=1
print("direct character and saturation bridges:",checks)

def boxsize(d,C):
 g=C-1
 v=(-1)**g+sum(2**i*comb(g+i,i) for i in range(2*d+1))
 den=2**(2*d+1)
 assert v%den==0
 return v//den

for d in range(1,8):
 for C in range(1,20):
  assert boxsize(d,C)==sum(
   comb(C-1-2*b+2*d-1,2*d-1) for b in range((C-1)//2+1))
print("profile count identities:",133)

for d in (7,8,10,12):
 C=3200*(d+1)**3+760*(d+1)**2
 print("profile box",d,"cutoff",C,"size",boxsize(d,C))

for r in range(1,13):
 d=2*r;a=4*r
 v=distance_value(d,[(1,1,a),(1,-1,a)])
 l=leading(d,2*a,0,0,0)
 assert v==comb(a,r)**2-comb(a,r-1)*comb(a,r)
 assert l==Q(a**(2*r),factorial(r)**2)
 ratio=Q(a-2*r+1,a-r+1)
 for i in range(r):ratio*=Q(a-i,a)**2
 assert v/l==ratio
print("linear-count Gaussian-margin family:",12)

def prior_test(d,cs):
 W=sum(n*c for n,e,c in cs);p=W-2*d
 if p<1:return True
 ep=(-1)**sum(c for n,e,c in cs if e<0)
 word=cs+[(p,ep,1)]
 nmax=max(n for n,e,c in word if c)
 delta0=(sum(n*c for n,e,c in word)-2*nmax)//2
 if delta0<=6 or nmax<=4:return True
 if all(e==1 for n,e,c in word if c):return True
 H=sum(c for n,e,c in word if n>=3)
 core=sum(n*c for n,e,c in word if n>=3)
 D=sum(n*c for n,e,c in word if n<3)
 if H==2 and core>=D:return True
 S=sum((n+1)**2*c for n,e,c in word if n>=5)
 J3=sum(c for n,e,c in word if n==3)
 J4=sum(c for n,e,c in word if n==4)
 V=sum(Q(c,2) if n==1 else c if (n,e)==(2,-1) else 0
       for n,e,c in word)
 bb=sum(c for n,e,c in word if (n,e)==(2,1))
 A=V+J3+J4;B=A+bb
 if S+16*J3+25*J4>=2**21*(nmax+1)**2:return True
 if S and B>=384*S and (A>=256 or B>=2**17):return True
 return False

assert 384*36-Q(5,2)==Q(27643,2)
assert 5*(2**21-2)//2==5242875
for d in (7,8,20,100):
 for m in (2*d+10,1000*d):
  assert not prior_test(d,[(1,1,m),(1,-1,m),(2,1,3)])
print("named-union filters: PASS")

def lpoly(d,t,b,k):
 H=hermite_polys(2*d,t-2*b);ans=[Q(0)]
 for j in range(d+1):
  ans=uniadd(ans,
   [Q(cat(j),factorial(d-j))*(t+k)**(d-j)*v for v in H[2*j]])
 return ans

def deval(poly,r,x):
 return sum(v*Q(factorial(i),factorial(i-r))*x**(i-r)
            for i,v in enumerate(poly) if i>=r)

checks=0
for d in range(1,7):
 for t,b,k in ((7,0,3),(0,5,2),(3,7,1),(0,0,11)):
  pp=[lpoly(j,t,b,k) for j in range(d+1)]
  lam=Q(2*(d+1),t+2*b+k)
  for x in (-10,0,11):
   ld=deval(pp[d],0,x)
   for h in range(d+1):
    for r in range(min(6,2*(d-h))+1):
     v=deval(pp[d-h],r,x)
     for aa in range(r+1):
      rhs2=4**r*(3*(d+1))**(2*aa)*lam**(2*h+r-aa)*ld**2
      assert (x**aa*v)**2<=rhs2
      checks+=1
print("differential-form exact checks:",checks)
print("ALL CHECKS PASS")