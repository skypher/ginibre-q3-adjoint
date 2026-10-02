import gzip, re, json, hashlib, random, os
from collections import defaultdict, Counter
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path

FM=Path('ginibre_q3/character_ring_iter/fm39')
def stamp(s): print(datetime.now(timezone.utc).isoformat(timespec='seconds'),s,flush=True)
def cg(a,b): return range(abs(a-b),a+b+1,2)
def table(w):
 d={(0,0):1}
 for z in sorted(w,key=lambda z:(-abs(z),z)):
  n=abs(z); eps=1 if z>0 else -1; q=defaultdict(int)
  for (r,s),v in d.items():
   for c in cg(r,n): q[c,s]+=v
   for c in cg(s,n): q[r,c]+=eps*v
  d={k:v for k,v in q.items() if v}
 return d
def cut(C):
 A=[];B=[];wa=wb=0
 for z in sorted(C,key=lambda z:(-abs(z),z)):
  if wa<=wb: A.append(z);wa+=abs(z)
  else: B.append(z);wb+=abs(z)
 if wa>wb:A,B=B,A
 return tuple(A),tuple(B)
def layers(L,i,j):
 A,B=cut(tuple(z for k,z in enumerate(L) if k not in (i,j)))
 x,y=table(A),table(B+(L[i],L[j])); o=defaultdict(int)
 for (r,s),v in x.items():
  if (r,s) in y:o[r+s]+=v*y[r,s]
 return dict(o)
def good(q): return all(v>=0 for v in q.values())
def pt(a,b): return tuple(sorted((a,b)))
def pos(L,p):
 for i in range(len(L)):
  for j in range(i+1,len(L)):
   if pt(L[i],L[j])==p:return i,j
 raise ValueError((L,p))
def candidates(L):
 out=[]
 def add(name,p):
  if p is not None and all(q[1]!=p for q in out):out.append((name,p))
 c=[]
 for i in range(len(L)):
  for j in range(i+1,len(L)):
   a,b=L[i],L[j]
   if (abs(a)-abs(b))%2==0:c.append(((abs(a)+abs(b),max(abs(a),abs(b)),pt(a,b),i,j),pt(a,b)))
 if c:add('top_same_parity',max(c)[1])
 neg=[i for i,z in enumerate(L) if z<0]
 if neg:
  i=max(neg,key=lambda i:(abs(L[i]),L[i],-i)); other=[j for j in range(len(L)) if j!=i]
  if other:add('largest_minus_other',pt(L[i],L[max(other,key=lambda j:(abs(L[j]),L[j],-j))]))
 ix=sorted(range(len(L)),key=lambda i:(abs(L[i]),L[i]),reverse=True)
 if len(ix)>=2:add('two_largest',pt(L[ix[0]],L[ix[1]]))
 ni=sorted(neg,key=lambda i:(abs(L[i]),L[i]),reverse=True)
 if len(ni)>=2:add('two_largest_minus',pt(L[ni[0]],L[ni[1]]))
 ix=sorted(range(len(L)),key=lambda i:(abs(L[i]),L[i]))
 if len(ix)>=2:add('minmax',pt(L[ix[0]],L[ix[-1]]))
 ps=sorted((i for i,z in enumerate(L) if z>0),key=lambda i:(abs(L[i]),L[i]))
 if len(ps)>=2:add('two_smallest_plus',pt(L[ps[0]],L[ps[1]]))
 return out
def choose(L):
 for name,p in candidates(L):
  q=layers(L,*pos(L,p))
  if good(q):return name,p,q
 return None,None,None
def read_rows(name):
 out=[]
 for s in gzip.decompress((FM/name).read_bytes()).decode().splitlines():
  if s.startswith('NOFLIP'):
   p=int(re.search(r'p=(-?\d+)',s).group(1))
   B=[int(x) for x in re.search(r'B=([-\d ]+?)\s+phi',s).group(1).split()]
   out.append(B+[p])
 return out
hard=[
 [1,1,-2]+[3]*6+[-4]*3+[5]*4+[6],
 [1,1,-2]+[3]*4+[-4]*3+[5]*5+[6],
 [1,1,1,-2,2]+[3]*5+[-4]*3+[4]*2+[5]*2+[8],
 [1,1,-2]+[3]*4+[-4]*3+[5]*4+[8],
 [-1]*6+[-2]*3+[-3]*3+[-4,-7],
 [-1,2,3,4,-5,6,7,8]]
rows40=read_rows('sec166_census_w40_noflip.log.gz')
rows52=read_rows('sec166_census_w52_noflip.log.gz')
base={tuple(sorted(x)):tuple(sorted(x,key=lambda z:(-abs(z),z))) for x in rows40+hard}
stamp(f'input W40={len(rows40)}, W52={len(rows52)}, hard=6, unique W40+hard={len(base)}')
which=Counter(); no_select=[]
for k,(key,L) in enumerate(sorted(base.items()),1):
 name,p,q=choose(L)
 if name is None:no_select.append(key)
 else:which[name]+=1
 if k%1000==0:stamp(f'candidate selector {k}/{len(base)}')
assert not no_select
stamp(f'finite candidate selector PASS; first-success counts={dict(which)}')
if os.environ.get('FM_STR9I_FULL','0')=='1':
 records=[];hist=Counter();ntypes=0
 for ix,(key,L) in enumerate(sorted(base.items()),1):
  works=set()
  for i,j in combinations(range(len(L)),2):
   if good(layers(L,i,j)):works.add(pt(L[i],L[j]))
  w=tuple(sorted(works));ntypes+=len(w);hist[len(w)]+=1
  records.append([key,w]);assert w,key
  if ix%500==0:stamp(f'all-pairs {ix}/{len(base)}')
 digest=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()
 stamp(f'pair-set PASS: {len(rows40)} census rows + 6 hard, {len(base)} unique multisets; working pair types={ntypes}; histogram={dict(sorted(hist.items()))}; sha256={digest}')
 for i,L in enumerate(hard):
  key=tuple(sorted(L));w=next(w for k,w in records if k==key)
  stamp(f'HARD[{i}] {key} => {w}')
 if os.environ.get('FM_STR9I_EMIT_ALL','0')=='1':
  for key,w in records:print('PAIRSET',key,w)
rng=random.Random(90210)
for L0 in rng.sample(rows52,300):
 assert choose(tuple(sorted(L0,key=lambda z:(-abs(z),z))))[0] is not None
stamp('W52 no-flip random sample 300 PASS')
rng=random.Random(314159)
for _ in range(300):
 n=rng.randint(3,8);L=[n,-n]+[(-1 if rng.randrange(2) else 1)*rng.randint(1,8) for __ in range(rng.randint(3,8))]
 if sum(z<0 for z in L)%2:
  j=next(j for j in range(2,len(L)) if L[j]!=0);L[j]=-L[j]
 assert choose(tuple(sorted(L,key=lambda z:(-abs(z),z))))[0] is not None
stamp('even-minus (+n,-n) random sample 300 PASS')
for k in range(6,19):
 w=k*(k+1)//2;p=next(p for p in range(k+1,k+4) if (w+p)%2==0)
 L=[-j for j in range(1,k+1)]+([p] if k%2==0 else [-p])
 assert choose(tuple(sorted(L,key=lambda z:(-abs(z),z))))[0]=='top_same_parity'
stamp('documented runs k=6..18, TopPair candidate PASS')
count=0
for n in range(1,9):
 for m in range(1,n+1):
  for q in range(1,n+1):
   for a in range(m,7):
    if (a-m)%2:continue
    L=[-n,-m,q]+[1]*a;i,j=pos(L,pt(-n,q))
    C=[z for t,z in enumerate(L) if t not in (i,j)];A,B=cut(C)
    assert sum(z<0 for z in A)==1 and sum(z<0 for z in B+(L[i],L[j]))==1
    assert sum(abs(z) for z in A)==sum(abs(z) for z in B)
    assert good(layers(L,i,j));count+=1
stamp(f'two-minus/+q/+1^a exact direct checks={count} PASS')
stamp('ALL CHECKS PASS')
