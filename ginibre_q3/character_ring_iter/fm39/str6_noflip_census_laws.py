from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import re
import os
import subprocess

PATH = "/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/fx6_52_noflip.log"
pat = re.compile(r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+) D=")
rows = 0
seen = set()
by_w = Counter()
by_l = Counter()
by_d = Counter()
by_mult = Counter()
by_ones = Counter()
by_odd = Counter()
by_changes = Counter()
LD = Counter()
max_l_d = {}
profiles = defaultdict(dict)
tp_min = None
hmin = None
hmax = None
monotone = alternating = one_sign = 0

for line in Path(PATH).open():
    m = pat.match(line)
    if not m:
        continue
    W, ps, bs, phi = int(m[1]), int(m[2]), m[3], int(m[4])
    B = tuple(map(int, bs.split()))
    p = ps
    assert sum(abs(z) for z in B) == W
    assert (W - abs(p)) % 2 == 0
    assert "D=ok" in line and re.search(r"TopPair\(.*\)=ok", line)
    word = B + (p,)
    count = Counter(abs(z) for z in word)
    sign_by_label = {}
    for z in word:
        n = abs(z)
        e = 1 if z > 0 else -1
        assert n not in sign_by_label or sign_by_label[n] == e
        sign_by_label[n] = e
    labels = sorted(sign_by_label)
    signs = tuple(sign_by_label[n] for n in labels)
    changes = sum(signs[i] != signs[i+1] for i in range(len(signs)-1))
    L = len(word)
    D = len(labels)
    mm = max(count.values())
    ones = count[1]
    odd = sum(count[n] for n in count if n % 2)
    delta = (W - abs(p)) // 2
    tp = re.search(r"TopPair\((-?\d+),(-?\d+)\)=ok", line)
    assert tp
    wt = abs(int(tp[1])) + abs(int(tp[2]))
    ratio = Fraction(wt, delta)
    h = sum((Fraction(1, n+1) for n in map(abs, word)), Fraction(0))
    rows += 1
    seen.add((W, p, B))
    by_w[W] += 1
    by_l[L] += 1
    by_d[D] += 1
    by_mult[mm] += 1
    by_ones[ones] += 1
    by_odd[odd] += 1
    by_changes[changes] += 1
    LD[L, D] += 1
    max_l_d[D] = max(max_l_d.get(D, 0), L)
    monotone += int(all(signs[i] <= signs[i+1] for i in range(D-1)) or
                    all(signs[i] >= signs[i+1] for i in range(D-1)))
    alternating += int(all(signs[i] != signs[i+1] for i in range(D-1)))
    one_sign += int(len(set(signs)) == 1)
    profile = tuple(sorted(count.items()))
    profiles[profile][signs] = (phi, W, p, B)
    if tp_min is None or ratio < tp_min[0]:
        tp_min = (ratio, W, p, B, wt, delta)
    if hmin is None or h < hmin[0]:
        hmin = (h, W, p, B)
    if hmax is None or h > hmax[0]:
        hmax = (h, W, p, B)

assert rows == 75532 and len(seen) == rows
assert all(L <= 2*D + 4 for L, D in LD)
assert max(L-D for L, D in LD) == 10
assert max(by_mult) == 4 and max(by_ones) == 2 and 2 not in by_odd
assert min(by_d) == 5 and max(by_d) == 10
assert max(by_l) == 16
assert len(profiles) == 26278
assert Counter(map(len, profiles.values())) == Counter({1:82, 2:16060, 3:10, 4:8941, 6:972, 8:213})
print("NOFLIP rows", rows, "weight counts", sorted(by_w.items()))
print("factor counts", sorted(by_l.items()))
print("distinct labels", sorted(by_d.items()), "max L by D", sorted(max_l_d.items()))
print("max multiplicity", sorted(by_mult.items()), "label-1 multiplicity", sorted(by_ones.items()))
print("odd-factor count", sorted(by_odd.items()), "sign transitions", sorted(by_changes.items()))
print("monotone sign sequences", monotone, "strict alternating", alternating, "one sign", one_sign)
print("absolute profiles", len(profiles), "no-flip signings per profile",
      sorted(Counter(map(len, profiles.values())).items()))
print("empirical L<=2D+4: max(L-2D) =", max(L-2*D for L, D in LD))
print("empirical L<=D+10: max(L-D) =", max(L-D for L, D in LD))
print("TopPair weight/delta minimum", tp_min)
print("sum 1/(n+1) range", hmin, hmax)

want = tuple((n, 1) for n in (1,2,3,4,5,6,7,9,10,11))
ex = profiles[want]
assert len(ex) == 8
print("same absolute labels, eight no-flip signings:", want)
for signs, (phi, W, p, B) in sorted(ex.items()):
    print("  signs", signs, "Phi", phi, "B", B, "sigma_p", p)


def cg(a,b):
    return range(abs(a-b),a+b+1,2)

def coeff2(factors,X,Y):
    rem=sum(abs(z) for z in factors)
    d={(0,0):1}
    for z in factors:
        n=abs(z);eps=1 if z>0 else -1;rem-=n
        out=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n):
                if abs(c-X)+abs(b-Y)<=rem:out[c,b]+=v
            for c in cg(b,n):
                if abs(a-X)+abs(c-Y)<=rem:out[a,c]+=eps*v
        d={k:v for k,v in out.items() if v}
    return d.get((X,Y),0)

profile_labels = tuple(n for n, _ in want)
for signs, (logged_phi, _, _, _) in ex.items():
    word = tuple(e*n for e, n in zip(signs, profile_labels))
    assert coeff2(word, 0, 0) == logged_phi

# Exact minimum over every admissible signing of this absolute-label profile.
base_labels = profile_labels[:-1]
profile_values = {}
for ss in product((-1, 1), repeat=len(base_labels)):
    psgn = -1 if sum(e < 0 for e in ss) % 2 else 1
    full_signs = ss + (psgn,)
    word = tuple(e*n for e, n in zip(full_signs, profile_labels))
    profile_values[full_signs] = coeff2(word, 0, 0)
profile_min = min(profile_values.values())
profile_minimizers = {sgn for sgn, val in profile_values.items() if val == profile_min}
assert len(profile_values) == 512 and profile_min == 50920
assert profile_minimizers == {
    (-1, 1, 1, 1, -1, 1, -1, -1, -1, -1),
    (1, 1, -1, 1, 1, 1, 1, 1, -1, 1),
}
assert profile_minimizers <= set(profiles[want])
assert len(profiles[want]) == 8
print("admissible signings", len(profile_values), "profile minimum", profile_min,
      "minimizing signings", sorted(profile_minimizers))

B=[1]*12+[2]*4+[3]*2
p=6
parent=coeff2(B+[p],0,0)
flipped=[-1,-1]+B[2:]+[p]
drop=parent-coeff2(flipped,0,0)
assert sum(B)==26 and (sum(B)-p)//2==10 and max(B)==3
assert parent==1188521234 and drop==1185770544 and drop%4==0
print("generic residual witness L=19 D=4, Phi=",parent,
      "flip(+1,+1) drop=",drop,"D=",drop//4)

Bodd = [3, 5] + [2]*9
odd_parent = coeff2(Bodd + [6], 0, 0)
odd_flip = coeff2([-3, -5] + [2]*9 + [6], 0, 0)
assert sum(Bodd) == 26 and (sum(Bodd)-6)//2 == 10
assert max(Bodd) == 5 and sum(n % 2 for n in Bodd+[6]) == 2
assert odd_parent == 639160 and odd_parent-odd_flip == 0
print("generic residual with two odd factors", "Phi=", odd_parent,
      "flip(+3,+5) drop=", odd_parent-odd_flip)


cpp = r'''#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using U = unsigned long long;
static int LIM;
static vector<int> c;
static U profiles=0, residuals=0;
static array<U,53> byW{}, byW_odd2{}, byW_bigmult{}, byW_many1{}, byW_long{};
static map<int,U> byL, byDistinct, byMaxMult, byOnes, byOdd, byChanges;
static U mono=0, alternating=0, oneSign=0, manyOddTwo=0, gt4=0, onesgt2=0, longL=0, L2Dfail=0, Ld10fail=0, joint=0;
static int sgn(int z){ return z>0?1:-1; }
static void emit(){
  ++profiles;
  if(profiles%10000000==0)cerr<<"profiles="<<profiles<<" residuals="<<residuals<<endl;
  int W=0,mx=0,fac=0,cores=0,neg=0,odd=0,ones=0,distinct=0,maxmult=0,negclasses=0;
  int prev=0,trans=0; bool alt=true,up=true,down=true,same=true;
  for(int n=1;n<=LIM;++n){
    int z=c[n],m=abs(z); if(!m)continue;
    W+=n*m;mx=n;fac+=m;distinct++;maxmult=max(maxmult,m);
    if(n>=3)cores+=m;
    if(z<0){neg+=m;negclasses++;}
    if(n==1)ones=m;
    if(n%2)odd+=m;
    int e=sgn(z);
    if(prev){if(e==prev)alt=false;else trans++;if(e<prev)up=false;if(e>prev)down=false;if(e!=prev)same=false;}
    prev=e;
  }
  if(cores<2)return;
  int sigma=(neg%2)?-1:1;
  for(int p=max(6,mx);p<=W;++p){
    if((W-p)%2)continue;
    int delta=(W-p)/2;
    if(delta<8||mx>delta)continue;
    if(p<=LIM&&c[p]*sigma<0)continue;
    ++residuals;++byW[W];
    bool newclass=p>mx;
    int L=fac+1;
    int d=distinct+(newclass?1:0);
    int mm=max(maxmult,newclass?1:abs(c[p])+1);
    int o=ones;
    int od=odd+(p%2);
    int tr=trans+(newclass&&prev!=sigma?1:0);
    bool al=alt&&(!newclass||prev!=sigma);
    bool monoUp=up&&(!newclass||prev<=sigma);
    bool monoDown=down&&(!newclass||prev>=sigma);
    bool one= same&&(!newclass||prev==sigma);
    byL[L]++;byDistinct[d]++;byMaxMult[mm]++;byOnes[o]++;byOdd[od]++;byChanges[tr]++;
    if(od==2){manyOddTwo++;byW_odd2[W]++;}
    if(mm>4)gt4++;if(o>2)onesgt2++;if(L>16)longL++;if(L>2*d+4)L2Dfail++;if(L>d+10)Ld10fail++;
    if(L<=16&&mm<=4&&o<=2&&od!=2)joint++;
    if(mm>4)byW_bigmult[W]++;
    if(o>2)byW_many1[W]++;
    if(L>16)byW_long[W]++;
    mono+=monoUp||monoDown;
    alternating+=al;
    oneSign+=one;
  }
}
static void gen(int n,int w,int fac,int cores,int neg,int odd,int ones,int maxmult,int distinct){
  if(n>LIM+1){emit();return;}
  int cap=(LIM-w)/n;
  for(int z=-cap;z<=cap;++z){
    int m=abs(z);
    c[n]=z;
    gen(n+1,w+n*m,fac+m,cores+(n>=3?m:0),neg+(z<0?m:0),
        odd+(n%2?m:0),ones+(n==1?m:0),max(maxmult,m),distinct+(m?1:0));
  }
  c[n]=0;
}
template<class M> static void dump(const char* name,const M& m){
  cout<<name;for(auto [k,v]:m)cout<<" "<<k<<":"<<v;cout<<"\n";
}
int main(int argc,char**argv){
  LIM=argc>1?atoi(argv[1]):52;if(LIM>52||LIM<1)return 2;
  c.assign(LIM+1,0);
  gen(1,0,0,0,0,0,0,0,0);
  U cases40=0,odd240=0,bigmult40=0,many140=0,long40=0;
  for(int w=1;w<=LIM;++w){
    if(w<=40){cases40+=byW[w];odd240+=byW_odd2[w];bigmult40+=byW_bigmult[w];many140+=byW_many1[w];long40+=byW_long[w];}
  }
  cout<<"LIM "<<LIM<<" profile_vectors "<<profiles<<" residual_cases "<<residuals
      <<" residual_W_le_40 "<<cases40<<" oddcount2_W_le_40 "<<odd240
      <<" maxmult_gt4_W_le_40 "<<bigmult40<<" ones_gt2_W_le_40 "<<many140
      <<" factors_gt16_W_le_40 "<<long40<<"\n";
  dump("factor_count",byL);dump("distinct_labels",byDistinct);dump("max_multiplicity",byMaxMult);
  dump("label1_count",byOnes);dump("odd_factor_count",byOdd);dump("sign_changes",byChanges);
  cout<<"sign_sequence_monotone "<<mono<<" strictly_alternating "<<alternating<<" one_sign "<<oneSign
      <<" odd_count_two "<<manyOddTwo<<"\n";
  cout<<"candidate_violations maxmult_gt4="<<gt4<<" ones_gt2="<<onesgt2
      <<" factors_gt16="<<longL<<" odd_count_eq2="<<manyOddTwo
      <<" L_gt_2D_plus_4="<<L2Dfail<<" L_gt_D_plus_10="<<Ld10fail
      <<" joint_filter_pass="<<joint<<endl;
}'''
fd = os.memfd_create("fmstr6_stats", 0)
os.fchmod(fd, 0o755)
env = dict(os.environ)
env["TMPDIR"] = "/dev/shm"
compiled = subprocess.run(
    ["g++", "-O3", "-std=c++17", "-x", "c++", "-o", f"/proc/self/fd/{fd}", "-"],
    input=cpp, text=True, env=env, pass_fds=(fd,), capture_output=True)
if compiled.returncode:
    print(compiled.stderr)
    raise SystemExit(compiled.returncode)
run = subprocess.run([f"/proc/self/fd/{fd}", "52"], env=env, pass_fds=(fd,),
                     text=True, capture_output=True)
print(run.stderr, end="")
print(run.stdout, end="")
assert run.returncode == 0
assert "residual_cases 516000611" in run.stdout
assert "maxmult_gt4=313755927" in run.stdout
assert "ones_gt2=327771730" in run.stdout
assert "factors_gt16=177560136" in run.stdout
assert "odd_count_eq2=15548502" in run.stdout
assert "L_gt_2D_plus_4=129135097" in run.stdout
assert "L_gt_D_plus_10=147573388" in run.stdout
assert "joint_filter_pass=120743304" in run.stdout
