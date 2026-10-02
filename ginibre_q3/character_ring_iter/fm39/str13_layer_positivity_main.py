# FM-STR13 (main agent): layer positivity of interior cuts.
# For a pair (u,v) of Lambda, cut C = Lambda - u - v by the interior rule (sort by decreasing |label|, give each factor
# to the lighter block, ties to A, swap so weight(A) <= weight(B)), put B' = B + {u,v}, and let
#   L_t = sum_{r+s=t} f_A(r,s) f_B'(r,s)   (the height-t layer; sum_t L_t = Phi(Lambda)).
# Checks on the W <= 40 no-flip census and on hard lists: (i) TopPair layer sign patterns; (ii) whether SOME pair has
# all layers >= 0; (iii) success counts of simple pair rules.
import re, gzip, pathlib, sys
from collections import defaultdict, Counter
def cg(a,b): return range(abs(a-b),a+b+1,2)
def table(word):
    d={(0,0):1}
    for z in sorted(word,key=abs,reverse=True):
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n): q[c,b]+=v
            for c in cg(b,n): q[a,c]+=e*v
        d={k:v for k,v in q.items() if v}
    return d
def interior_cut(C):
    A=[];B=[];wa=wb=0
    for z in sorted(C,key=lambda z:(-abs(z),z)):
        if wa<=wb: A.append(z); wa+=abs(z)
        else: B.append(z); wb+=abs(z)
    if wa>wb: A,B=B,A
    return A,B
def toppair(B):
    best=None
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if (abs(B[i])-abs(B[j]))%2: continue
            key=(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])))
            if best is None or key>best[0]: best=(key,i,j)
    return best[1],best[2]
def layers(L,i,j):
    C=[z for q,z in enumerate(L) if q not in (i,j)]; A,Bb=interior_cut(C)
    fa=table(A); fb=table(Bb+[L[i],L[j]]); lay=defaultdict(int)
    for k,v in fa.items():
        if k in fb: lay[k[0]+k[1]]+=v*fb[k]
    return lay
def allpos(lay): return all(v>=0 for v in lay.values())
LOG=pathlib.Path(__file__).resolve().parent/"sec166_census_w40_noflip.log.gz"
rows=[]
for line in gzip.decompress(LOG.read_bytes()).decode().splitlines():
    if not line.startswith("NOFLIP"): continue
    p=int(re.search(r"p=(-?\d+)",line).group(1)); Bk=[int(x) for x in re.search(r"B=([-\d ]+?)\s+phi",line).group(1).split()]
    rows.append(Bk+[p])
hard=[[1,1,-2,3,3,3,3,3,3,-4,-4,-4,5,5,5,5,6],[1,1,-2,3,3,3,3,3,-4,-4,-4,5,5,5,5,5,6],[1,1,1,-2,2,3,3,3,3,3,-4,-4,-4,4,4,5,5,8],
      [1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8],[-1]*6+[-2]*3+[-3]*3+[-4]+[-7],[-1,2,3,4,-5,6,7,8]]
pat=Counter(); no_pair=[]
for L in rows+hard:
    i,j=toppair(L[:-1]); lay=layers(L,i,j)
    s=''.join('+' if lay[t]>0 else '-' for t in sorted(lay) if lay[t])
    pat[''.join(c for k,c in enumerate(s) if k==0 or c!=s[k-1])]+=1
    if not allpos(lay):
        if not any(allpos(layers(L,a,b)) for a in range(len(L)) for b in range(a+1,len(L))): no_pair.append(L)
print("lists",len(rows)+len(hard),"TopPair layer sign patterns",dict(pat))
print("lists with NO pair having all layers >= 0:",len(no_pair),no_pair[:3])
assert not no_pair
print("FM-STR13 LAYER POSITIVITY (SOME PAIR) PASS")
