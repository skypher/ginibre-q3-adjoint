# FM-STR11 main-agent check: on every no-flip list of the W <= 40 census, no split with both halves of size >= MINSIZE
# (default 2) has f_A f_B >= 0 in every channel.  Usage: python3 -u str11_noflip_split_kill.py [MINSIZE] [SAMPLE]
import itertools, sys, time, re, random, gzip, pathlib
from functools import lru_cache
from collections import defaultdict
def cg(a,b): return range(abs(a-b),a+b+1,2)
@lru_cache(None)
def table(word):
    d={(0,0):1}
    for z in word:
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n): q[c,b]+=v
            for c in cg(b,n): q[a,c]+=e*v
        d={k:v for k,v in q.items() if v}
    return d
def compatible(A,B):
    fa=table(A); fb=table(B)
    return all(v*fb.get(k,0)>=0 for k,v in fa.items())
def splits(word, minsize):
    cls=sorted(set(word)); mult=[word.count(c) for c in cls]
    for take in itertools.product(*[range(m+1) for m in mult]):
        A=tuple(sorted(sum(([c]*t for c,t in zip(cls,take)),[])))
        B=tuple(sorted(sum(([c]*(m-t) for c,m,t in zip(cls,mult,take)),[])))
        if len(A)<minsize or len(B)<minsize or A>B: continue
        yield A,B
rows=[]
LOG=pathlib.Path(__file__).resolve().parent/'sec166_census_w40_noflip.log.gz'
MINSIZE=int(sys.argv[1]) if len(sys.argv)>1 else 2
for line in gzip.decompress(LOG.read_bytes()).decode().splitlines():
    if not line.startswith("NOFLIP"): continue
    p=int(re.search(r"p=(-?\d+)",line).group(1)); B=[int(x) for x in re.search(r"B=([-\d ]+?)\s+phi",line).group(1).split()]
    rows.append(tuple(sorted(B+[p])))
rows=sorted(set(rows)); random.seed(1); 
if len(sys.argv)>2: rows=random.sample(rows,int(sys.argv[2]))
print('minimum half size',MINSIZE)
print(time.strftime('%H:%M:%S'),"distinct no-flip lists",len(rows),flush=True)
stat=defaultdict(int); bad=[]; t=time.time()
for i,w in enumerate(rows):
    best=None
    for A,B in splits(w,MINSIZE):
        if compatible(A,B): best=(A,B); break
    if best is None: bad.append(w)
    else: stat[min(len(best[0]),len(best[1]))]+=1
    if time.time()-t>30: print(time.strftime('%H:%M:%S'),i+1,"/",len(rows),"no compatible split:",len(bad),flush=True); t=time.time()
print(time.strftime('%H:%M:%S'),"done",len(rows),"failures",len(bad),"first-found min half sizes",dict(stat))
print("examples",bad[:10])
assert len(bad)==len(rows), 'a no-flip list has a compatible split'
print('FM-STR11 NO-FLIP SPLIT CHECK PASS')
