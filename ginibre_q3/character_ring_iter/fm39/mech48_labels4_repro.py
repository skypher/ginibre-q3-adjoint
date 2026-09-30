import argparse, os, subprocess, random
from functools import lru_cache
from fractions import Fraction as Q
from math import factorial, comb

parser = argparse.ArgumentParser(
    description="FM-MECH48: exact sector-4 certificate, memory-only C++/GMP")
parser.add_argument("--threads", type=int, default=4)
args = parser.parse_args()

rho = Q(7, 9)
delta = Q(2, 9)
tau = Q(19, 27)
J = Q(901, 4050)

@lru_cache(None)
def gb(N, d):
    K = N*delta+d+1
    return max(
        (K+b)*(K+b+1) *
        (Q(1, 4*(b+1)*(b+2))+J*tau**b)/delta**2
        for b in range(5))

@lru_cache(None)
def eb(N, b):
    low = Q(factorial(N), 2**(N+1))
    for j in range(1, N+2):
        low /= b+j
    high = (rho**(N+1)-Q(2, 5)**(N+1))/Q(N+1)*tau**b
    return (low+high)/delta**2

def sup(l, m, n):
    if m:
        return 3*Q(7, 10)**l*Q(5, 8)**(m-1)*Q(2, 3)**n
    if l >= 2:
        return 4*Q(7, 10)**(l-2)*Q(2, 3)**n
    if l:
        return 2*Q(2, 3)**n
    return Q(2, 3)**n

H37 = Q(5629, 200)*Q(7, 10)**35*Q(749, 9)*Q(758, 9)
assert H37 < 1
assert Q(7, 10)*Q(38, 37)**2 < 1
print("H>=37 bound:", H37)

states = []
count = 0
maxN = maxB = 0
for H in range(1, 37):
    for l in range(H+1):
        for m in range(H-l+1):
            n = H-l-m
            if m+n == 0:
                continue
            R = max(1, (l+1)//2+m)
            S = sup(l, m, n)
            raw = Q(2)**l*Q(3)**m*Q(2, 3)**n
            d = l+m+2*n
            for N in range(max(2, R), 100):
                if S*rho**(N-R)*gb(N, d) < 1:
                    break
                K = N*delta+d+1

                def good(b):
                    return raw*(K+b)*(K+b+1)*eb(N, b) < 1

                hi = 4
                while not good(hi):
                    hi *= 2
                lo = 3
                while hi > lo+1:
                    mid = (lo+hi)//2
                    if good(mid):
                        hi = mid
                    else:
                        lo = mid
                B = hi
                assert good(B)
                ways = sum(
                    max(0, N-max(1, (al+m+1)//2)
                        -(l-al+m+1)//2+1)
                    for al in range(l+1))
                if not ways:
                    continue
                states.append((N, l, m, n, B))
                count += ways*B
                maxN = max(maxN, N)
                maxB = max(maxB, B)
            else:
                raise AssertionError((l, m, n))
    if H % 6 == 0:
        print("H", H, "certified states", len(states), "entries", count)
print("Box states", len(states), "entries", count,
      "Nmax", maxN, "Bmax", maxB)

CPP = r'''
#include <gmpxx.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
using Z=mpz_class;
using Table=std::vector<std::vector<Z>>;
constexpr int H=36, QH=37;
struct Plan {
 int nmax=-1,jmax=0,tmax=0;
 std::array<int,QH> mm,ll,dd;
 Plan(){mm.fill(-1);ll.fill(-1);dd.fill(-1);}
};
size_t ix(int T,int m,int r){
 return size_t(m)*(T+1)-size_t(m)*(m-1)/2+r;
}
std::vector<std::vector<Z>> moments;
const Z& M(int m,int r,int b){
 return moments[m+r+b][ix(m+r+b,m,r)];
}
void exact_ui(Z& x,unsigned long d){
 assert(mpz_divisible_ui_p(x.get_mpz_t(),d));
 mpz_divexact_ui(x.get_mpz_t(),x.get_mpz_t(),d);
}
void make_moments(int TMAX){
 std::vector<Z> cat(TMAX+1);cat[0]=1;
 for(int j=1;j<=TMAX;++j){
  cat[j]=cat[j-1]*(4*j-2);exact_ui(cat[j],j+1);
 }
 moments.resize(TMAX+1);
 for(int T=0;T<=TMAX;++T){
  moments[T].resize((T+1)*(T+2)/2);
  for(int m=0;m<=T/2;++m){
   int r=T-m,deg=2*T,d=2*(m-r);
   Z c0=1,c1=d,ans=cat[T],tmp;
   for(int j=1;j<deg;++j){
    tmp=d*c1-(deg-j+1)*c0;exact_ui(tmp,j+1);
    c0=c1;c1=tmp;
    if((j+1)%2==0)
     ans+=c1*cat[(j+1)/2]*cat[T-(j+1)/2];
   }
   moments[T][ix(T,m,r)]=ans;
   moments[T][ix(T,r,m)]=ans;
  }
  for(int b=1;b<=T;++b)for(int m=0;m<=T-b;++m){
   int r=T-b-m;
   Z v=moments[T][ix(T,m+1,r)]+moments[T][ix(T,m,r+1)];
   exact_ui(v,2);
   v-=2*moments[T-1][ix(T-1,m,r)];
   moments[T][ix(T,m,r)]=std::move(v);
  }
  if(T%25==0){
   std::printf("Moment degree %d/%d\n",T,TMAX);
   std::fflush(stdout);
  }
 }
}
Table mixed(int m,int r,int J,int D){
 Table t(J+1);
 for(int j=0;j<=J;++j)t[j].resize(D-j+1);
 for(int b=0;b<=D;++b)t[0][b]=M(m,r,b);
 for(int j=1;j<=J;++j)for(int b=0;b<=D-j;++b){
  Z v=8*(m-r)*t[j-1][b];
  if(j>=2)
   v+=4*(j-1)*t[j-2][b+1]+8*(j-1)*t[j-2][b];
  if(b)v+=12*b*t[j][b-1];
  exact_ui(v,2*(m+r+b+j+2));
  t[j][b]=std::move(v);
 }
 return t;
}
Z shifted(int m,int r,int j,int b){
 Z ans=0,c,term;
 for(int h=0;h<=j;++h){
  mpz_bin_uiui(c.get_mpz_t(),j,h);
  term=c*M(m+h,r+j-h,b);
  if((j-h)%2)ans-=term;else ans+=term;
 }
 assert(mpz_divisible_2exp_p(ans.get_mpz_t(),2*j));
 mpz_tdiv_q_2exp(ans.get_mpz_t(),ans.get_mpz_t(),2*j);
 return ans;
}
int main(int argc,char**argv){
 if(argc>1&&(!std::strcmp(argv[1],"-h")||
              !std::strcmp(argv[1],"--help"))){
  std::puts("FM-MECH48 exact box verifier; "
            "reads certified box and direct checks on stdin");
  return 0;
 }
 int NMAX,nrec,threads;
 unsigned long long expected;
 std::cin>>NMAX>>nrec>>expected>>threads;
 omp_set_num_threads(threads);
 auto id=[](int N,int l,int m,int n){
  return ((size_t(N)*QH+l)*QH+m)*QH+n;
 };
 std::vector<unsigned char> bound(
     size_t(NMAX+1)*QH*QH*QH,0);
 std::vector<Plan> plans(NMAX+1);
 for(int z=0;z<nrec;++z){
  int N,l,m,n,B;std::cin>>N>>l>>m>>n>>B;
  assert(B<256);
  bound[id(N,l,m,n)]=B;
  auto&p=plans[N];p.nmax=std::max(p.nmax,n);
  p.mm[n]=std::max(p.mm[n],m);
  p.ll[n]=std::max(p.ll[n],l);
  p.dd[n]=std::max(p.dd[n],B-1+l+m);
  p.jmax=std::max(p.jmax,l+2*n);
  p.tmax=std::max(p.tmax,B-1+l+m+2*n);
 }
 int nq;std::cin>>nq;
 using Key=std::array<int,7>;
 std::map<Key,Z> queries;
 for(int z=0;z<nq;++z){
  Key key;for(int&i:key)std::cin>>i;
  std::string value;std::cin>>value;
  queries.emplace(key,Z(value));
 }
 assert(int(queries.size())==nq);
 int TMAX=0;
 std::vector<std::array<int,2>> jobs;
 for(int N=2;N<=NMAX;++N)if(plans[N].nmax>=0){
  TMAX=std::max(TMAX,N+plans[N].tmax);
  for(int a=0;a<N;++a)jobs.push_back({N,a});
 }
 std::printf(
   "Box entries %llu; backgrounds %zu; moment cutoff %d\n",
   expected,jobs.size(),TMAX);
 std::fflush(stdout);
 make_moments(TMAX);
 std::vector<std::vector<std::vector<std::pair<int,long>>>>
     coeff(QH);
 for(int l=0;l<=H;++l){
  coeff[l].resize(l+1);
  for(int a=0;a<=l;++a){
   std::vector<long> row(l+1);row[0]=1;
   if(l)row[1]=2*a-l;
   for(int j=1;j<l;++j){
    long v=(2*a-l)*row[j]-(l-j+1)*row[j-1];
    assert(v%(j+1)==0);row[j+1]=v/(j+1);
   }
   for(int j=0;j<=l;++j)
    if(row[j])coeff[l][a].push_back({j,row[j]});
  }
 }
 unsigned long long total=0,zeros=0,bridges=0,seen=0;
 int done=0;Z least=0;
 double started=omp_get_wtime();
 #pragma omp parallel for schedule(dynamic)
 for(size_t job=0;job<jobs.size();++job){
  int N=jobs[job][0],a=jobs[job][1],r=N-a,A=2*a,E=2*r;
  const Plan&p=plans[N];
  Table base=mixed(a,r,p.jmax,p.tmax);
  unsigned long long local=0,zlocal=0,blocal=0,qlocal=0;
  Z minlocal=0,value;
  for(int j: {1,p.jmax})if(j<=p.jmax){
   for(int b: {0,p.tmax-j}){
    assert(base[j][b]==shifted(a,r,j,b));++blocal;
   }
  }
  for(int n=0;n<=p.nmax;++n){
   if(p.mm[n]>=0){
    int J=p.ll[n],D=p.dd[n];
    Table work(J+1);
    for(int j=0;j<=J;++j)
     work[j]=std::vector<Z>(
         base[j].begin(),base[j].begin()+D-j+1);
    for(int m=0;m<=p.mm[n];++m){
     for(int l=0;l<=J;++l){
      int B=bound[id(N,l,m,n)];
      if(!B)continue;
      int lo=std::max(0,l+m-A),hi=std::min(l,E-m);
      for(int al=lo;al<=hi;++al)for(int b=0;b<B;++b){
       value=0;
       for(auto [j,c]:coeff[l][al]){
        const Z&v=work[j][b+l-j];
        if(c>0)
         mpz_addmul_ui(value.get_mpz_t(),v.get_mpz_t(),c);
        else
         mpz_submul_ui(value.get_mpz_t(),v.get_mpz_t(),-c);
       }
       if(value<0){
        #pragma omp critical
        {
         std::cerr<<"NEGATIVE "<<A<<" "<<E<<" "<<b<<" "
          <<al<<" "<<l-al<<" "<<m<<" "<<n<<" "<<value<<"\n";
         std::cerr.flush();
        }
        std::exit(2);
       }
       ++local;
       if(value==0)++zlocal;
       else if(minlocal==0||value<minlocal)minlocal=value;
       if(N<=4&&b<=2&&l+m+n<=4){
        Key key{A,E,b,al,l-al,m,n};
        auto it=queries.find(key);
        if(it!=queries.end()){
         assert(value==it->second);++qlocal;
        }
       }
      }
     }
     if(m<p.mm[n])for(int j=0;j<=J;++j)
      for(int b=0;b<D-m-j;++b)
       mpz_sub(work[j][b].get_mpz_t(),
               work[j][b+1].get_mpz_t(),
               work[j][b].get_mpz_t());
    }
   }
   if(n<p.nmax)for(int j=0;j<=p.jmax-2*n-2;++j)
    for(int b=0;b<=p.tmax-2*n-2-j;++b){
     mpz_add(base[j][b].get_mpz_t(),
             base[j][b+2].get_mpz_t(),
             base[j][b+1].get_mpz_t());
     mpz_submul_ui(base[j][b].get_mpz_t(),
                  base[j+2][b].get_mpz_t(),2);
    }
  }
  #pragma omp critical
  {
   total+=local;zeros+=zlocal;bridges+=blocal;seen+=qlocal;
   if(minlocal>0&&(least==0||minlocal<least))least=minlocal;
   ++done;
   if(done%16==0||done==int(jobs.size())){
    std::printf("Backgrounds %d/%zu; words %llu; seconds %.1f\n",
                done,jobs.size(),total,
                omp_get_wtime()-started);
    std::fflush(stdout);
   }
  }
 }
 assert(total==expected);assert(seen==queries.size());
 std::cout<<"EXACT BOX "<<total<<" zeros "<<zeros
          <<" least "<<least
          <<"\nSHIFT BRIDGES "<<bridges
          <<"\nDIRECT BRIDGES "<<seen<<"\nPASS\n";
}
'''

def add(*terms):
    out = {}
    for c, poly in terms:
        for key, v in poly.items():
            out[key] = out.get(key, 0)+c*v
    return {key: v for key, v in out.items() if v}

def mul(p, q):
    out = {}
    for (i, j), v in p.items():
        for (k, l), w in q.items():
            key = i+k, j+l
            out[key] = out.get(key, 0)+v*w
    return {key: v for key, v in out.items() if v}

def U(n, axis):
    return {
        ((n-2*j, 0) if axis == 0 else (0, n-2*j)):
        (-1)**j*comb(n-j, j)
        for j in range(n//2+1)}

ONE = {(0, 0): 1}
S = {(1, 0): 1, (0, 1): 1}
D = {(1, 0): 1, (0, 1): -1}
Z = {(2, 0): 1, (0, 2): 1, (0, 0): -2}
P = {(1, 1): 1}

def hs(k):
    out = {}
    for j in range(k+1):
        out = add((1, out), (1, mul(U(j, 0), U(k-j, 1))))
    return out

H2 = hs(2)
H3 = hs(3)
S3 = add((1, U(3, 0)), (1, U(3, 1)))
S4 = add((1, U(4, 0)), (1, U(4, 1)))
assert H2 == add((1, Z), (1, P))
assert H3 == mul(S, add((1, Z), (-1, ONE)))
assert S3 == mul(S, add((1, Z), (-1, P)))
assert S4 == add((1, mul(Z, Z)), (1, Z), (-2, mul(P, P)))

GEN = {"s": S, "d": D, "z": Z, "h2": H2,
       "h3": H3, "S3": S3, "S4": S4}

@lru_cache(None)
def power(key, n):
    return ONE if n == 0 else mul(power(key, n-1), GEN[key])

@lru_cache(None)
def cat(n):
    return comb(2*n, n)//(n+1)

def direct(key):
    A, E, b, al, ga, m, n = key
    suffix = A-ga-m
    out = ONE
    for g, ex in (
        ("d", E), ("s", suffix), ("z", b), ("h2", al),
        ("h3", m), ("S3", ga), ("S4", n)
    ):
        out = mul(out, power(g, ex))
    return sum(
        v*cat(i//2)*cat(j//2)
        for (i, j), v in out.items() if i % 2 == j % 2 == 0)

pool = []
for N, l, m, n, B in states:
    if N > 4 or l+m+n > 4:
        continue
    for al in range(l+1):
        ga = l-al
        lo = (ga+m+1)//2
        hi = N-max(1, (al+m+1)//2)
        for aa in range(lo, hi+1):
            for b in range(min(B, 3)):
                pool.append((2*aa, 2*(N-aa), b, al, ga, m, n))

keys = random.Random(48).sample(pool, 200)
checks = [(key, direct(key)) for key in keys]
print("Direct Catalan expectations prepared:", len(checks))

packet = f"{maxN} {len(states)} {count} {args.threads}\n"
packet += "".join(" ".join(map(str, s))+"\n" for s in states)
packet += str(len(checks))+"\n"
packet += "".join(
    " ".join(map(str, key))+" "+str(v)+"\n"
    for key, v in checks)

def locate(name):
    return subprocess.check_output(
        ["g++", "-print-file-name="+name], text=True).strip()

obj = os.memfd_create("mech48-object", 0)
exe = os.memfd_create("mech48-executable", 0)
subprocess.run(
    ["g++", "-O3", "-std=c++17", "-fno-pie", "-fopenmp",
     "-pipe", "-x", "c++", "-", "-c", "-o",
     f"/proc/self/fd/{obj}"],
    input=CPP, text=True, pass_fds=(obj,), check=True)
subprocess.run(
    ["ld", "-L"+os.path.dirname(locate("crtbegin.o")),
     "-o", f"/proc/self/fd/{exe}",
     "-dynamic-linker", "/lib64/ld-linux-x86-64.so.2",
     locate("crt1.o"), locate("crti.o"), locate("crtbegin.o"),
     f"/proc/self/fd/{obj}", locate("libstdc++.so"),
     "-lgmpxx", "-lgmp", locate("libgomp.so"), "-lm",
     locate("libgcc_s.so"), "-lc",
     locate("crtend.o"), locate("crtn.o")],
    pass_fds=(obj, exe), check=True)
subprocess.run(
    [f"/proc/self/fd/{exe}"], input=packet, text=True,
    pass_fds=(exe,), check=True)
os.close(obj)
os.close(exe)
