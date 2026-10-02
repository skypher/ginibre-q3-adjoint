# FM-STR15 (main agent): Lancaster-coupling positivity D(rho) = sum_r f(r,r) rho^r >= 0 on [-1,1],
# certified exactly by Bernstein coefficients with subdivision.  Small box (labels <= 4, length <= 7) and the W <= 40 census.
# D(rho) = sum_r f(r,r) rho^r >= 0 on [-1,1]?  Exact Bernstein test with subdivision; W<=40 census + small box.
import re, gzip, itertools, sys, time
from fractions import Fraction as Fr
from math import comb
import pathlib
_here=pathlib.Path(__file__).resolve().parent
src=open(_here/"str13_layer_positivity_main.py").read(); src=src[:src.index("LOG=pathlib")]
ns={}; exec(src,ns); table=ns['table']
def bern_nonneg(c, lo=Fr(-1), hi=Fr(1), depth=0):
    # c: power-basis coefficients in rho; test p>=0 on [lo,hi] via Bernstein coefficients after affine map
    d=len(c)-1
    # p(lo + (hi-lo) t) coefficients in t
    a=[Fr(0)]*(d+1)
    for k,ck in enumerate(c):
        if ck==0: continue
        # (lo + w t)^k
        w=hi-lo
        for j in range(k+1): a[j]+=ck*comb(k,j)*lo**(k-j)*w**j
    b=[sum(a[i]*Fr(comb(j,i),comb(d,i)) for i in range(j+1)) for j in range(d+1)]
    if min(b)>=0: return True
    # endpoint negative -> real failure
    if b[0]<0 or b[-1]<0: return False
    if depth>12: return None
    mid=(lo+hi)/2
    r1=bern_nonneg(c,lo,mid,depth+1); r2=bern_nonneg(c,mid,hi,depth+1)
    if r1 is False or r2 is False: return False
    if r1 is None or r2 is None: return None
    return True
def D(word):
    f=table(word); R=max([r for r,s in f]+[0]); return [f.get((r,r),0) for r in range(R+1)]
LOG=_here/"sec166_census_w40_noflip.log.gz"
rows=[]
for line in gzip.decompress(open(LOG,'rb').read()).decode().splitlines():
    if not line.startswith("NOFLIP"): continue
    p=int(re.search(r"p=(-?\d+)",line).group(1)); B=[int(x) for x in re.search(r"B=([-\d ]+?)\s+phi",line).group(1).split()]
    rows.append(B+[p])
labs=[s*n for n in range(1,5) for s in (1,-1)]; small=[]
for N in range(2,8):
    for w in itertools.combinations_with_replacement(sorted(labs),N):
        if sum(1 for z in w if z<0)%2 or sum(map(abs,w))%2: continue
        small.append(list(w))
for name,lists in (("small box labels<=4 len<=7",small),("W<=40 no-flip census",rows)):
    t0=time.time(); bad=[]; unk=0
    for w in lists:
        r=bern_nonneg(D(w))
        if r is False: bad.append(w)
        elif r is None: unk+=1
    print(time.strftime('%H:%M:%S'),name,len(lists),"failures",len(bad),"undecided",unk,bad[:4],f"{time.time()-t0:.0f}s",flush=True)
