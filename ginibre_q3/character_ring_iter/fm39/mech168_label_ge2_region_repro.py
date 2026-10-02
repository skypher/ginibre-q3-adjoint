import argparse,re,time
from pathlib import Path
from math import comb
from fractions import Fraction as Q
from functools import lru_cache
from collections import Counter
from itertools import combinations
from random import Random
ap=argparse.ArgumentParser(description="FM-MECH168 exact verifier; read-only, Python integers.")
ap.add_argument("--skip-census",action="store_true")
ap.add_argument("--census",type=Path,default=Path(
 "/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/"
 "scratchpad/flipdesc/fx3b_48.log"))
args=ap.parse_args();started=time.monotonic()

def cg(a,b):return range(abs(a-b),a+b+1,2)
@lru_cache(None)
def fusion(ns):
 v=[1]
 for n in ns:
  w=[0]*(len(v)+n+2)
  for i,c in enumerate(v):
   if c:w[abs(i-n)]+=c;w[i+n+2]-=c
  for i in range(2,len(w)):w[i]+=w[i-2]
  v= w[:-2]
 return tuple(v)
def inv(ns):return fusion(tuple(sorted(ns)))[0]
def feasible(ns):
 return not ns or (len(ns)>=2 and sum(ns)%2==0 and 2*max(ns)<=sum(ns))

# Profiles are scaled by 6, so capacities use integer dot products.
@lru_cache(None)
def profile(ns,negative=False,old=False):
 ns=tuple(sorted(ns));N=len(ns)
 if not ns:return (6,)
 if not feasible(ns):return ()
 if N==2:return (6,)*(ns[0]+1)
 if N==3:
  S=sum(ns)//2
  return tuple(6*max(0,1+min(ns[0],2*r,S-r,S-ns[-1]+r))
               for r in range(S+1))
 v=[6]*(ns[-1]+1);m=ns[-3]
 for r in range(m+1):
  v[r]=max(v[r],6*(r+1) if 2*r<=m else 6) if old else max(v[r],3*min(2*r+2,m+2))
 if negative and N==4 and any(c%2 for c in Counter(ns).values()):
  assert min(ns)>=2
  for r,c in enumerate((6,16,20,18)):v[r]=max(v[r],c)
 return tuple(v)
@lru_cache(None)
def capacity(A,B,negative,old):
 return sum(a*b for a,b in zip(profile(A,negative,old),profile(B,negative,old)))
def cost(counts):return sum((Q(36*c,k) for k,c in counts.items()),Q(0))
def top(B):
 pairs=(ij for ij in combinations(range(len(B)),2) if (B[ij[0]]-B[ij[1]])%2==0)
 return max(pairs,key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                                max(abs(B[ij[0]]),abs(B[ij[1]]))))
def subsets(ns):
 return [tuple(ns[i] for i in range(len(ns)) if s>>i&1)
         for s in range(1<<len(ns))]
def masks(s):
 t=s
 while True:
  yield t
  if not t:return
  t=(t-1)&s

def certificate(B,p,mode=2):
 L=len(B)+1;ns=tuple(map(abs,B))+(abs(p),)
 F=(1<<L)-1;bg=F>>1;sg=sum(1<<i for i,z in enumerate(B) if z<0)
 A=subsets(ns);i,j=top(B);child=F^(1<<i)^(1<<j)
 neg=Counter();pos=Counter();old=(mode==0)
 for s in range(1,bg+1):
  if (s&sg).bit_count()%2 and feasible(A[s]) and feasible(A[F^s]):
   neg[capacity(A[s],A[F^s],True,old)]+=1
 for s in masks(child&bg):
  if s and not (s&sg).bit_count()%2 and feasible(A[s]) and feasible(A[child^s]):
   pos[capacity(A[s],A[child^s],False,old)]+=1
 gain=Q(min(ns[i],ns[j])+1)
 if mode==2 and feasible(A[child]):
  h=profile(A[child])
  gain=Q(sum(h[r] for r in range(abs(ns[i]-ns[j])//2,(ns[i]+ns[j])//2+1)),6)
 rho=1-cost(neg)-(1+cost(pos))/gain if feasible(A[child]) else 1-cost(neg)
 return rho,gain,A,child,i,j

# Independent profile checks on genuine fusion rows.
rng=Random(168);checks=0
for _ in range(1200):
 ns=tuple(sorted(rng.randrange(2,15) for _ in range(rng.randrange(2,9))))
 if not feasible(ns):continue
 f=fusion(ns);negative=any(c%2 for c in Counter(ns).values())
 for r,h in enumerate(profile(ns,negative)):
  assert 6*(f[2*r] if 2*r<len(f) else 0)>=h*f[0];checks+=1
assert fusion((3,3,6))[2]==2*fusion((3,3,6))[0]
print("profile checks",checks,"PASS; sharper 9/4 candidate fails at (3,3,6)")

def capacities(m):
 q=m//2
 if m%2:return Q((q+1)*(3*q+5),2),Q((q+1)*(2*q+3)*(8*q+13),12)
 return Q((q+1)*(3*q+2),2),Q((q+1)*(8*q*q+13*q+6),6)
def nu(L):return 2**(L-2)-(L if L%2==0 else 2*L-2)
def bad(L):return L-2+3*comb(L-2,2)
def margin(L,m):
 aa,bb=capacities(m);J=L-2;K=2**(J-1)-1-J-comb(J,2)
 q=min(nu(L),bad(L));h=min(K,bad(J))
 beta=Q(q)/aa+Q(nu(L)-q)/bb
 P=Q(1,6)+Q(comb(J,2)+h)/aa+Q(K-h)/bb
 return 1-beta-(1+P)/3
expected=[19,23,27,32,38,45,54,65]
for L,M0 in zip(range(9,17),expected):
 M=4
 while margin(L,M)<0:M+=1
 assert M==M0
 print("two-small threshold",L,M)
for L in range(9,65):
 lo,hi=0,1
 while hi**3<7*2**L:hi*=2
 while lo+1<hi:
  v=(lo+hi)//2
  if v**3>=7*2**L:hi=v
  else:lo=v
 assert margin(L,max(4*L,hi)-1)>=0

caller_count=0
for L,M in ((9,19),(12,32),(16,65)):
 for small in ((2,2),(2,3),(3,3)):
  B=list(small)+[M+2*(i%2) for i in range(L-3)]
  B=[-n if n%3 else n for n in B]
  p=max(map(abs,B))+1
  if (sum(map(abs,B))-p)%2:p+=1
  rho,gain,A,child,i,j=certificate(tuple(B),p)
  F=(1<<L)-1;bg=F>>1;sg=sum(1<<i for i,z in enumerate(B) if z<0)
  ms=[inv(a) for a in A]
  g=sum((-1 if (s&sg).bit_count()%2 else 1)*ms[s]*ms[F^s] for s in range(bg+1))
  hv=sum((-1 if (s&sg).bit_count()%2 else 1)*ms[s]*ms[child^s] for s in masks(child&bg))
  assert rho>=margin(L,M)>=0 and ms[F]>0
  assert Q(g-hv)>=rho*ms[F]
  caller_count+=1
print("two-small caller checks",caller_count,"PASS")

# The flip witnesses use a two-variable Clebsch-Gordan evaluator.
@lru_cache(None)
def table(B):
 v={(0,0):1}
 for z in B:
  n=abs(z);eps=1 if z>0 else -1;w={}
  for (i,j),c in v.items():
   for k in cg(i,n):w[k,j]=w.get((k,j),0)+c
   for k in cg(j,n):w[i,k]=w.get((i,k),0)+eps*c
  v={ij:c for ij,c in w.items() if c}
 return v
def flip(full,i,j):
 C=tuple(z for k,z in enumerate(full) if k not in (i,j))
 return (1 if full[j]>0 else -1)*table(C).get((abs(full[i]),abs(full[j])),0)
def U(n):return {n-2*j:(-1)**j*comb(n-j,j) for j in range(n//2+1)}
def mul(P,Q0):
 R={}
 for (i,j),c in P.items():
  for (k,l),d in Q0.items():R[i+k,j+l]=R.get((i+k,j+l),0)+c*d
 return {ij:c for ij,c in R.items() if c}
def factor(n,e):
 P={}
 for i,c in U(n).items():
  P[i,0]=P.get((i,0),0)+c;P[0,i]=P.get((0,i),0)+e*c
 return {ij:c for ij,c in P.items() if c}
def word(ns):
 P={(0,0):1}
 for z in ns:P=mul(P,factor(abs(z),1 if z>0 else -1))
 return P

for q,h,p in ((2,4,8),(3,5,10)):
 A=128;cores=(-q,-q,-(q+4),-(q+4),p)
 B=(h,)*A+cores[:-1];full=B+(p,);i,j=top(B)
 assert {i,j}=={A+2,A+3}
 polys={"parent":word(cores),"child":word((-q,-q,p))}
 exps={"parent":A,"child":A}
 for v in (-q,-(q+4),p,h):
  C=list(cores);C.remove(-q)
  power=A-1 if v==h else A
  if v!=h:C.remove(v)
  xy={(i,j):a*b for i,a in U(q).items() for j,b in U(abs(v)).items()}
  P=mul(word(C),xy)
  polys[str(v)]={ij:(1 if v>0 else -1)*c for ij,c in P.items()}
  exps[str(v)]=power
 polys["other"]=mul(word(cores),{(i,j):a*b for i,a in U(h).items() for j,b in U(h).items()})
 exps["other"]=A-2
 J=max(max(i,j) for P in polys.values() for i,j in P)
 ballot=[{r:comb(j,(j-r)//2)-(comb(j,(j-r)//2-1) if j-r>=2 else 0)
          for r in range(j%2,j+1,2)} for j in range(J+1)]
 mom=[];v=[1]
 for a in range(A+1):
  mom.append([sum(c*v[r] for r,c in row.items() if r<len(v)) for row in ballot])
  w=[0]*(len(v)+h+2)
  for r,c in enumerate(v):
   if c:w[abs(r-h)]+=c;w[r+h+2]-=c
  for r in range(2,len(w)):w[r]+=w[r-2]
  v=w[:-2]
 @lru_cache(None)
 def joint(a,i,j):
  return sum(comb(a,z)*mom[z][i]*mom[a-z][j] for z in range(a+1))
 def value(P,a):return sum(c*joint(a,i,j) for (i,j),c in P.items())
 for a in range(5):
  assert value(polys["parent"],a)==table((h,)*a+cores).get((0,0),0)
 vals={s:value(P,exps[s]) for s,P in polys.items()}
 assert vals["parent"]%2==vals["child"]%2==0
 g0,h0=vals["parent"]//2,vals["child"]//2
 assert h0<2*g0 and 3*g0<2*h0
 assert max(vals[str(v)] for v in (-q,-(q+4),p,h))<0<vals["other"]
 W=sum(map(abs,B));delta=(W-p)//2
 assert W-p==2*delta and p>=max(map(abs,B)) and max(map(abs,B))<=delta
 print("small-label obstruction",q,"A=128; 1/2 < parent/child < 2/3;",
       "all small flips < 0 < repeated-background flip; W",W,"delta",delta)

if not args.skip_census:
 rows=[]
 for line in args.census.read_text().splitlines():
  m=re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+)",line)
  if not m:continue
  B=tuple(map(int,m[3].split()))
  if len(B)+1>=9 and min(map(abs,B)) in (2,3):
   rows.append((B,int(m[2]),int(m[4])))
 assert len(rows)==832
 counts=[Counter() for _ in range(3)];miss=[]
 for num,(B,sp,ph) in enumerate(rows,1):
  p=abs(sp);L=len(B)+1;F=(1<<L)-1;bg=F>>1
  sg=sum(1<<i for i,z in enumerate(B) if z<0)
  assert (-1 if sg.bit_count()%2 else 1)*p==sp
  ns=tuple(map(abs,B))+(p,);state={}
  for z in B+(sp,):
   assert abs(z) not in state or state[abs(z)]==(z>0)
   state[abs(z)]=(z>0)
  for mode in range(3):
   rho,gain,A,child,i,j=certificate(B,p,mode)
   if rho>=0:counts[mode][(L,min(ns))]+=1
  ms=[inv(a) for a in A]
  w=[(-1 if (s&sg).bit_count()%2 else 1)*ms[s]*ms[F^s] for s in range(bg+1)]
  g=sum(w);assert 2*g==ph
  hv=sum((-1 if (s&sg).bit_count()%2 else 1)*ms[s]*ms[child^s] for s in masks(child&bg))
  assert ms[F]>=gain*ms[child]
  assert Q(g-hv)>=rho*ms[F]
  nz=[(s,v) for s,v in enumerate(w) if v]
  for u,v in combinations(range(L),2):
   assert sum(t for s,t in nz if ((s>>u)^(s>>v))&1)<0
  T,X,Y,Z=[sum(v for s,v in nz if 2*((s>>i)&1)+((s>>j)&1)==t) for t in range(4)]
  assert 2*g==2*T+(X+Y)+(Y+Z)+(X+Z)
  assert T-hv-abs(X)-abs(Y)-abs(Z)>=0
  if rho<0:miss.append((B,sp))
  if num%200==0:print("census callers",num,"/832")
 assert [sum(c.values()) for c in counts]==[422,786,798]
 assert counts[2]==Counter({(9,2):726,(9,3):10,(10,2):62})
 for mode,c in enumerate(counts):print("census region",mode,sorted(c.items()),"total",sum(c.values()))
 templates={
 (2,3,4,5,5,6,6,7,8):{10,18,20},
 (2,3,4,5,5,6,7,7,8):{9,19,21},
 (2,3,4,5,5,6,7,7,9):{10,12,14,16,18,20,22},
 (2,3,4,5,5,6,7,8,8):{10,18,20,22}}
 assert len(miss)==34
 for B,sp in miss:
  ns=tuple(map(abs,B));assert abs(sp) in templates[ns]
  assert all(z<0 for z in B+(sp,)) or all((z<0)==(abs(z)%2==0) for z in B+(sp,))
 print("complement: 34 lists, four templates and odd-label reflections")
print("FM-MECH168 PASS; seconds",round(time.monotonic()-started,3))
