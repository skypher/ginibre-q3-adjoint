import argparse,itertools,random
from collections import Counter,defaultdict
from functools import lru_cache
from math import comb

ap=argparse.ArgumentParser(description='Exact FM-CHK99 independent censuses; run from repository root.')
ap.add_argument('--seven-max',type=int,default=16)
ap.add_argument('--word-max',type=int,default=12)
ap.add_argument('--eight-max',type=int,default=9)
ap.add_argument('--progress-step',type=int,default=5000)
args=ap.parse_args()

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def inv(ns):
    ns=tuple(sorted(ns)); k=len(ns)
    if k==0:return 1
    total=sum(ns)
    if total&1 or 2*max(ns)>total:return 0
    if k==1:return 0
    if k==2:return int(ns[0]==ns[1])
    if k==3:return 1
    if k==4:
        lo=max(abs(ns[0]-ns[1]),abs(ns[2]-ns[3]))
        hi=min(ns[0]+ns[1],ns[2]+ns[3])
        return max(0,(hi-lo)//2+1)
    out=0
    for mask in range(1<<k):
        t=total//2-sum(ns[i]+1 for i in range(k) if mask>>i&1)+k-2
        if t>=k-2:
            out+=(-1)**mask.bit_count()*comb(t,k-2)
    assert out>=0
    return out

def inv_fusion(ns):
    row={0:1}
    for n in ns:
        nxt=defaultdict(int)
        for x,v in row.items():
            for y in cg(x,n):nxt[y]+=v
        row=nxt
    return row.get(0,0)

def sign_masks(ns,pair_free=False):
    groups=[];i=0
    while i<len(ns):
        j=i+1
        while j<len(ns) and ns[j]==ns[i]:j+=1
        groups.append((i,j-i));i=j
    choices=[(0,m) if pair_free else tuple(range(m+1)) for _,m in groups]
    for ks in itertools.product(*choices):
        mask=0;neg=0
        for (start,m),k in zip(groups,ks):
            neg+=k
            for i in range(start,start+k):mask|=1<<i
        yield mask,neg

def phi_subset(word):
    ns=tuple(sorted(abs(x) for x in word)); f=len(ns)
    negmask=sum(1<<i for i,x in enumerate(word) if x<0)
    if negmask.bit_count()&1:return 0
    if sum(ns)&1:return 0
    total=inv(ns)
    for I in itertools.combinations(range(f),2):
        m=(1<<I[0])|(1<<I[1])
        sgn=-1 if (negmask&m).bit_count()&1 else 1
        comp=tuple(ns[k] for k in range(f) if not(m>>k&1))
        total+=sgn*inv((ns[I[0]],ns[I[1]]))*inv(comp)
    for r in (3,4):
        if r>f//2:continue
        for I in itertools.combinations(range(f),r):
            if 2*r==f and 0 not in I:continue
            m=sum(1<<k for k in I)
            sgn=-1 if (negmask&m).bit_count()&1 else 1
            A=tuple(ns[k] for k in I)
            B=tuple(ns[k] for k in range(f) if not(m>>k&1))
            total+=sgn*inv(A)*inv(B)
    return 2*total

def phi_direct(word):
    row={(0,0):1}
    for z in word:
        n=abs(z); eps=1 if z>0 else -1
        nxt=defaultdict(int)
        for (a,b),v in row.items():
            for x in cg(a,n):nxt[x,b]+=v
            for y in cg(b,n):nxt[a,y]+=eps*v
        row=nxt
    return row.get((0,0),0)

rng=random.Random(991273)
for _ in range(300):
    ns=tuple(rng.randint(1,12) for _ in range(rng.randint(1,8)))
    assert inv(ns)==inv_fusion(ns),(ns,inv(ns),inv_fusion(ns))
print('multiplicity formula vs CG DP: 300 PASS',flush=True)

# Theorem 2 budget: pair-free seven-factor words, labels 2..L.
L=args.seven_max
assert L>=2
profiles=checked=0; min_pay=None; min_count=None; min_N=None; min_slack_word=None
for ns in itertools.combinations_with_replacement(range(2,L+1),7):
    profiles+=1
    N=inv(ns)
    P=sum(inv(tuple(ns[k] for k in range(7) if k not in (i,j)))
          for i,j in itertools.combinations(range(7),2) if ns[i]==ns[j])
    tdata=[]
    for I in itertools.combinations(range(7),3):
        mask=sum(1<<i for i in I)
        A=tuple(ns[i] for i in I)
        B=tuple(ns[i] for i in range(7) if not(mask>>i&1))
        tdata.append((mask,inv(A)*inv(B)))
    for negmask,negcount in sign_masks(ns,True):
        if negcount&1:continue
        vals=[w for m,w in tdata if (negmask&m).bit_count()&1]
        d=max(vals,default=0);T=sum(vals)
        s1=N+P-20*d;s2=20*d-T
        if s1<0 or s2<0:
            word=tuple(-n if negmask>>i&1 else n for i,n in enumerate(ns))
            exact=phi_subset(word);direct=phi_direct(word)
            print('COUNTEREXAMPLE',word,'N,P,d,T,slacks',N,P,d,T,s1,s2,
                  'Phi',exact,'direct',direct,flush=True)
            raise SystemExit(1)
        checked+=1
        if min_pay is None or s1<min_pay:
            min_pay=s1
            min_slack_word=tuple(-n if negmask>>i&1 else n for i,n in enumerate(ns))
        min_count=s2 if min_count is None else min(min_count,s2)
        min_N=N if min_N is None else min(min_N,N)
    if profiles%args.progress_step==0:
        print('seven-budget progress',profiles,'/',comb(L+6,7),
              'words',checked,flush=True)
print('seven-budget L=',L,'profiles=',profiles,'pair-free signings=',checked,
      'min(N+P-20d)=',min_pay,'at',min_slack_word,
      'min(20d-T)=',min_count,'min N=',min_N,'PASS',flush=True)

# Theorem 3: all signed seven-factor multisets, labels 1..H.
H=args.word_max
all_words=evald=forced_oddminus=forced_oddweight=0
min_phi=None;min_word=None
for ns in itertools.combinations_with_replacement(range(1,H+1),7):
    total_signings=1
    for _,m in Counter(ns).items():total_signings*=m+1
    all_words+=total_signings
    if sum(ns)&1:
        forced_oddweight+=total_signings
        continue
    for negmask,negcount in sign_masks(ns,False):
        if negcount&1:
            forced_oddminus+=1
            continue
        word=tuple(-n if negmask>>i&1 else n for i,n in enumerate(ns))
        v=phi_subset(word)
        if v<0:
            d=phi_direct(word)
            print('THEOREM3 COUNTEREXAMPLE',word,'Phi',v,'direct',d,flush=True)
            raise SystemExit(1)
        evald+=1
        if min_phi is None or v<min_phi:min_phi=v;min_word=word
    if all_words%args.progress_step<total_signings:
        print('seven-all progress words=',all_words,'/',comb(2*H+6,7),
              'evaluated=',evald,flush=True)
print('seven-all H=',H,'signed multisets=',all_words,'evaluated=',evald,
      'forced odd total weight=',forced_oddweight,
      'forced odd minus parity=',forced_oddminus,
      'min Phi=',min_phi,'at',min_word,'PASS',flush=True)

# Theorem 4: all signed eight-factor words in the stated sectors.
J=args.eight_max
sector_counts=Counter();sector_eval=Counter();sector_min={};sector_minword={}
for ns in itertools.combinations_with_replacement(range(2,J+1),8):
    if min(ns)>=3:
        sector='all labels >=3'
    elif min(ns)>=2 and any(n&1 for n in ns):
        sector='min=2, some odd'
    else:
        continue
    sector_counts[sector]+=1
    for negmask,negcount in sign_masks(ns,False):
        if negcount&1 or sum(ns)&1:continue
        word=tuple(-n if negmask>>i&1 else n for i,n in enumerate(ns))
        v=phi_subset(word)
        if v<0:
            d=phi_direct(word)
            print('THEOREM4 COUNTEREXAMPLE',word,'Phi',v,'direct',d,flush=True)
            raise SystemExit(1)
        sector_eval[sector]+=1
        if sector not in sector_min or v<sector_min[sector]:
            sector_min[sector]=v;sector_minword[sector]=word
print('eight sectors J=',J,'profiles=',dict(sector_counts),
      'parity-admissible signed cases=',dict(sector_eval),
      'minima=',dict(sector_min),'at=',dict(sector_minword),'PASS',flush=True)

# Independent two-variable checks on deterministic samples.
for f,Hh in ((7,args.word_max),(8,args.eight_max)):
    nsrange=range(1,Hh+1) if f==7 else range(2,Hh+1)
    samples=[]
    for ns in itertools.combinations_with_replacement(nsrange,f):
        if f==8 and not (min(ns)>=3 or
                         (min(ns)>=2 and any(n&1 for n in ns))):continue
        if sum(ns)&1:continue
        for negmask,negcount in sign_masks(ns,False):
            if negcount&1:continue
            samples.append(tuple(-n if negmask>>i&1 else n
                                 for i,n in enumerate(ns)))
            if len(samples)>=80:break
        if len(samples)>=80:break
    for word in samples:
        assert phi_subset(word)==phi_direct(word),(word,phi_subset(word),phi_direct(word))
    print('two-variable evaluator cross-check f=',f,
          'cases=',len(samples),'PASS',flush=True)
