import argparse,hashlib,itertools,pathlib,re,subprocess
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import comb

ap=argparse.ArgumentParser(description="Exact FM-MECH153 checks.")
ap.add_argument("--max-label",type=int,default=7)
ap.add_argument("--receipts",action="store_true")
ap.add_argument("--census-log")
args=ap.parse_args()

def cg(a,b): return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def fusion(ns):
    d={0:1}
    for n in ns:
        e=Counter()
        for a,v in d.items():
            for b in cg(a,n): e[b]+=v
        d=dict(e)
    return d

@lru_cache(None)
def inv(ns):
    k=len(ns)
    if not k:return 1
    s=sum(ns)
    if s%2 or 2*max(ns)>s:return 0
    if k==1:return 0
    if k==2:return int(ns[0]==ns[1])
    if k==3:return 1
    if k==4:
        a,b,c,d=ns
        return max(0,(min(a+b,c+d)-max(abs(a-b),abs(c-d)))//2+1)
    out=0
    for mask in range(1<<k):
        t=s//2-sum(ns[i]+1 for i in range(k) if mask>>i&1)+k-2
        if t>=k-2:out+=(-1)**mask.bit_count()*comb(t,k-2)
    return out

def cuts(ns):
    f=len(ns);out=[]
    for k in range(2,f//2+1):
        for ii in itertools.combinations(range(f),k):
            m=sum(1<<i for i in ii)
            if 2*k==f and not m&1:continue
            a=tuple(ns[i] for i in range(f) if m>>i&1)
            b=tuple(ns[i] for i in range(f) if not m>>i&1)
            out.append((m,k,inv(a)*inv(b)))
    return out

def signs(ns):
    classes=sorted(set(ns))
    for bits in range(1<<len(classes)):
        m=sum(1<<i for i,n in enumerate(ns)
              if bits>>classes.index(n)&1)
        if m.bit_count()%2==0:yield m

def bands(I,J):
    return tuple(sum(t in cg(a,b) for a in I for b in J)
                 for t in (0,2,4,6))

intervals=[tuple(range(a,b+1,2)) for a in range(25)
           for b in range(a,25,2)]
bandchecks=0
for I in intervals:
    for J in intervals:
        K=set(I)&set(J);d=len(K)
        if not d or min(K)==0 or I==J:continue
        if len(I)>=3 and len(J)>=3:
            low=((1,3,6,7) if d==1 else
                 (2,6,7,6) if d==2 else
                 (d,3*d-1,5*d-5,7*d-12))
            assert all(x>=y for x,y in zip(bands(I,J),low))
            bandchecks+=1
        if I[0]%2==J[0]%2==0 and len(I)%2==len(J)%2==0:
            low=((1,3,3,1) if d==1 else
                 (2,5,7,5) if d==2 else
                 (3,9,12,12) if d==3 else
                 (d,3*d-1,5*d-5,7*d-12))
            assert all(x>=y for x,y in zip(bands(I,J),low))
            bandchecks+=1
print("interval endpoint checks",bandchecks,"PASS")

E=(F(1),F(5,2),F(3),F(2))
M=(F(1),F(5,2),F(5,2),F(3,2))
O=(F(1),F(5,2),F(3),F(1))
Os=[(F(1),F(3),F(3),F(1)),
    (F(1),F(5,2),F(7,2),F(5,2)),
    (F(1),F(3),F(4),F(4)),
    (F(1),F(11,4),F(15,4),F(4))]
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
assert min(dot(a,b) for a in Os for b in Os)==20
assert dot(E,E)>=20 and dot(E,M)>=10
assert dot(E,O)>=12 and dot(M,M)>=12 and dot(M,O)>=10

quadchecks=0
for ns in itertools.combinations_with_replacement(range(1,13),4):
    d=inv(ns)
    if not d:continue
    f=fusion(ns);o=sum(n%2 for n in ns)
    classes=sorted(set(ns))
    canneg=any(sum(ns.count(n) for i,n in enumerate(classes)
                   if mask>>i&1)%2 for mask in range(1<<len(classes)))
    if not canneg:continue
    v=E if o==0 else M if o==2 else O
    assert all(f.get(t,0)>=d*x for t,x in zip((0,2,4,6),v))
    if o==4:
        v=Os[min(d,4)-1]
        assert all(f.get(t,0)>=d*x for t,x in zip((0,2,4,6),v))
    quadchecks+=1
print("mixed-sign quartet profiles",quadchecks,"PASS")

for o in range(0,9,2):
    ceiling={0:20,2:10,4:12,6:10,8:20}[o]
    for m in range(0,9,2):
        for r in range(max(0,m-8+o),min(m,o)+1):
            odd=set(range(o))
            minus=set(range(r))|set(range(o,o+m-r))
            v=sum(len(set(S)&odd)%2==0 and
                  len(set(S)&minus)%2==1
                  for S in itertools.combinations(range(8),4))//2
            assert v<=ceiling
print("eight-factor parity ceilings PASS")

checked7=checked8=0
H=args.max_label
for f in (7,8):
    for ns in itertools.combinations_with_replacement(
            range(2 if f==7 else 1,H+1),f):
        if sum(ns)%2:continue
        N=inv(ns);cs=cuts(ns)
        for sg in signs(ns):
            P=0;d3=d4=0;t4=0
            for mask,k,w in cs:
                neg=(mask&sg).bit_count()%2
                if k==2:
                    assert not neg or w==0
                    P+=w
                if k==3 and neg:d3=max(d3,w)
                if k==4 and neg:d4=max(d4,w);t4+=w
            if f==7:
                assert N+P>=20*d3
                checked7+=1
            else:
                assert N>=t4
                if all(n%2 for n in ns):assert N>=20*d4
                checked8+=1
print("seven/eight signed checks",checked7,checked8,"PASS")

def statistics(word):
    word=tuple(sorted(word,key=abs));ns=tuple(map(abs,word))
    sg=sum(1<<i for i,n in enumerate(word) if n<0)
    out=Counter(N=inv(ns),P=0,n3=0,p3=0,n4=0,p4=0)
    for mask,k,v in cuts(ns):
        neg=(mask&sg).bit_count()%2
        if k==2:out["P"]+=(-1)**neg*v
        else:out[("n" if neg else "p")+str(k)]+=v
    out["phi"]=2*(out["N"]+out["P"]+out["p3"]+out["p4"]
                  -out["n3"]-out["n4"])
    return out

for w in ((1,1,-2,-2,3,3,3,5),(-1,2,3,4,-5,6,7,8)):
    a=statistics(w);diff=[]
    for i,j in itertools.combinations(range(8),2):
        v=list(w);v[i]*=-1;v[j]*=-1
        diff.append(a["phi"]-statistics(tuple(v))["phi"])
    print("budget test",w,dict(a),"flip difference range",min(diff),max(diff))
    if w==(-1,2,3,4,-5,6,7,8):
        assert a==dict(N=553,P=0,n3=111,p3=60,n4=66,p4=42,phi=956)
        assert (min(diff),max(diff))==(-708,-8)
        assert a["P"]+a["p3"]+a["p4"]-a["n3"]==-9
assert inv((1,1,1,1,1,1,3,3))==20
assert statistics((-1,-1,-1,-1,-1,-1,3,3))["n4"]==20

if args.receipts:
    root=pathlib.Path("ginibre_q3")
    expected={
      "verify_su2_seven_shallow_z3.cpp":"5c4f65cd45ca6d5237f43d0641b97cf288edc6a2f251405f605e9b2325ed86aa",
      "verify_su2_seven_shallow_rank_two_cells_z3.cpp":"a55ea6da6a01e2c07a4ed0e72745493f72808fd0ad792d1c290a8f3ed1f2d0e8",
      "collect_su2_seven_shallow_exact.cpp":"e173390c1a611b46f76da6fd807765c6aade844ff74dba1c07733f49facba372",
      "verify_su2_d12_deep_minus_z3.cpp":"32a05b998980ac167f7595effbe9e93bfc31858ae1c7dd42a170d5a6cfd60a52"}
    for n,h in expected.items():
        assert hashlib.sha256((root/"character_ring_iter"/n).read_bytes()).hexdigest()==h
    rr=root/"certificates/su2_seven_shallow_exact_receipts"
    lines=(rr/"SHA256SUMS").read_text().splitlines()
    for line in lines:
        h,n=line.split(maxsplit=1)
        assert hashlib.sha256(pathlib.Path(n.strip().lstrip("*")).read_bytes()).hexdigest()==h
    subprocess.run([str(root/"character_ring_iter/collect_su2_seven_shallow_exact"),
        str(root/"character_ring_iter/verify_su2_seven_shallow_z3"),
        str(rr/"rank_one_aggregate.log"),str(rr)],check=True)
    log=(root/"certificates/su2_d12_deep_minus_z3.log").read_text()
    assert "queries=768" in log and "counterexamples=UNSAT result=PASS" in log
    print("source identities and",len(lines),"receipt hashes PASS")

if args.census_log:
    lengths=Counter();t3zero=0
    for line in pathlib.Path(args.census_log).read_text().splitlines():
        m=re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+)",line)
        if not m:continue
        w=tuple(map(int,m[3].split()))+(int(m[2]),)
        lengths[len(w)]+=1
        if len(w)==8:
            a=statistics(w)
            assert a["phi"]==int(m[4])
            t3zero+=a["n3"]==0
    assert lengths=={6:499,7:1597,8:1944,9:1066,10:270,11:52,12:2}
    assert t3zero==0
    print("no-flip census",dict(sorted(lengths.items())),
          "eight-factor zero negative 3|5:",t3zero)
print("FM-MECH153 PASS")
