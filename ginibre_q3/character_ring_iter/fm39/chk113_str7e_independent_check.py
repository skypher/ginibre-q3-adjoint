import argparse, gzip, os, re, shlex, subprocess, random
from pathlib import Path

ap = argparse.ArgumentParser(description="Fresh exact FM-STR7e checker; anonymous in-memory C++ build.")
ap.add_argument("--threads", type=int, default=24)
ap.add_argument("--census", default="/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39/sec166_census_w40_noflip.log.gz")
args = ap.parse_args()
assert args.threads >= 1

def stamp(s):
    from datetime import datetime, timezone
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), s, flush=True)

pat = re.compile(r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+) D=")
rows = []
for line in gzip.open(args.census, "rt"):
    m = pat.match(line)
    if not m:
        continue
    W = int(m.group(1))
    p = int(m.group(2))
    B = tuple(map(int, m.group(3).split()))
    phi = int(m.group(4))
    assert sum(abs(x) for x in B) == W
    rows.append((B + (p,), phi))
assert len(rows) == 5430
stamp("loaded gz census rows=5430")

jobs = []
def add_job(group, word, expected=-1, pair=(-1,-1), note=""):
    jobs.append((group, tuple(word), expected, pair[0], pair[1], note))

for word, phi in rows:
    add_job(0, word, phi)

runs = []
for k in range(1,21):
    W = k*(k+1)//2
    for p in range(k+1,k+11):
        delta_num = W-p
        if p < max(6,k) or delta_num % 2:
            continue
        delta = delta_num//2
        if delta < max(8,k) or k > delta:
            continue
        if sum(i >= 3 for i in range(1,k+1)) < 2:
            continue
        sig = -1 if k % 2 else 1
        w = tuple(-i for i in range(1,k+1)) + (sig*p,)
        assert sum(abs(x) for x in w[:-1]) == W
        runs.append(w)
        add_job(1, w)
assert len(runs) == 68
stamp("generated admissible run words=68")

# Three explicit matching receipts: the requested sign minimum and two census rows.
minimum = (-1,2,3,4,-5,6,7,8)
add_job(2, minimum, 956, (4,6), "specified-sign-minimum")
for ix in (0,1):
    w, phi = rows[ix]
    pos = [j for j,x in enumerate(w) if x == -5]
    assert len(pos) >= 2
    add_job(2, w, phi, (pos[0],pos[1]), "census-row-"+str(ix+1))

# Deterministic random signed words for the two identities.
rr = random.Random(113112)
random_jobs = 250
for q in range(random_jobs):
    L = rr.randint(3,10)
    labels = [rr.randint(1,9) for _ in range(L)]
    signs = [1 if rr.randrange(2) else -1 for _ in range(L)]
    if sum(s<0 for s in signs)%2: signs[0] *= -1
    w = tuple(s*n for s,n in zip(signs,labels))
    i,j = sorted(rr.sample(range(L),2))
    add_job(3, w, -1, (i,j), "random-"+str(q))

CPP = r"""
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;

static void stamp(const std::string& s) {
    std::time_t t=std::time(nullptr); std::tm u; gmtime_r(&t,&u);
    std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
struct Poly {
    int top;
    std::vector<Z> v;
    explicit Poly(int d=0): top(d), v((d+1)*(d+1)) {}
    Z get(int x,int y) const {
        if(x<0||y<0||x+y>top) return 0;
        return v[x*(top+1)+y];
    }
    Z& at(int x,int y) { return v[x*(top+1)+y]; }
};
static Poly axis_product(const Poly& f,int n,int axis) {
    Poly g(f.top+n);
    for(int fixed=0;fixed<=f.top;++fixed) {
        int edge=f.top-fixed;
        std::vector<Z> p0(edge+2),p1(edge+2);
        for(int u=0;u<=edge;++u) {
            Z z=axis==0?f.get(u,fixed):f.get(fixed,u);
            p0[u+1]=p0[u]; p1[u+1]=p1[u];
            (u&1?p1:p0)[u+1]+=z;
        }
        for(int t=0;t+fixed<=g.top;++t) {
            int lo=std::abs(t-n), hi=std::min(edge,t+n);
            int parity=(t+n)&1;
            if((lo&1)!=parity) ++lo;
            if(lo>hi) continue;
            Z z=(parity?p1:p0)[hi+1]-(parity?p1:p0)[lo];
            if(axis==0) g.at(t,fixed)+=z;
            else g.at(fixed,t)+=z;
        }
    }
    return g;
}
static Poly add(const Poly& a,const Poly& b,int sign=1) {
    assert(a.top==b.top); Poly c(a.top);
    for(size_t i=0;i<c.v.size();++i)c.v[i]=a.v[i]+sign*b.v[i];
    return c;
}
static Poly signed_factor(const Poly& f,int signed_label) {
    int n=std::abs(signed_label),sg=signed_label>0?1:-1;
    Poly x=axis_product(f,n,0), y=axis_product(f,n,1);
    return add(x,y,sg);
}
static int mass(const std::vector<int>& w) {
    int s=0; for(int x:w)s+=std::abs(x); return s;
}
static Poly product(std::vector<int> w) {
    std::stable_sort(w.begin(),w.end(),[](int a,int b){return std::abs(a)>std::abs(b);});
    Poly p(0); p.at(0,0)=1;
    for(int z:w)p=signed_factor(p,z);
    return p;
}
static Poly one_term(Poly p,const std::vector<std::pair<int,int>>& ops) {
    for(auto [n,axis]:ops)p=axis_product(p,n,axis);
    return p;
}
static void add_scaled(Poly& dst,const Poly& src,int c) {
    assert(dst.top==src.top);
    for(size_t i=0;i<dst.v.size();++i)dst.v[i]+=c*src.v[i];
}
static Poly pair_pure(const Poly& base,int a,int b,int eu,int ev) {
    Poly out(base.top+a+b);
    Poly q=one_term(base,{{a,0},{b,0}}); add_scaled(out,q,1);
    q=one_term(base,{{a,1},{b,1}}); add_scaled(out,q,eu*ev);
    return out;
}
static Poly pair_mixed(const Poly& base,int a,int b,int eu,int ev) {
    Poly out(base.top+a+b);
    Poly q=one_term(base,{{a,0},{b,1}}); add_scaled(out,q,ev);
    q=one_term(base,{{a,1},{b,0}}); add_scaled(out,q,eu);
    return out;
}
static std::pair<std::vector<int>,std::vector<int>> interior_split(std::vector<int> c) {
    std::stable_sort(c.begin(),c.end(),[](int x,int y){return std::abs(x)>std::abs(y);});
    std::vector<int> a,b; int wa=0,wb=0;
    for(int x:c) {
        if(wa<=wb){a.push_back(x);wa+=std::abs(x);}
        else {b.push_back(x);wb+=std::abs(x);}
    }
    if(wa>wb) std::swap(a,b);
    return {a,b};
}
static std::vector<int> remainder(const std::vector<int>& w,int i,int j) {
    std::vector<int> c; c.reserve(w.size()-2);
    for(int k=0;k<(int)w.size();++k)if(k!=i&&k!=j)c.push_back(w[k]);
    return c;
}
static Z dot_height(const Poly& a,const Poly& b,int h) {
    Z s=0;
    for(int r=0;r<=h;++r)s+=a.get(r,h-r)*b.get(r,h-r);
    return s;
}
struct PairData {
    std::vector<Z> P,M;
    Z phi=0,D=0;
    int first_budget_bad=-1,first_plain_bad=-1;
    int prefix_count=0;
};
static PairData evaluate_pair(const std::vector<int>& w,int i,int j,const Z& phi) {
    int a=std::abs(w[i]),b=std::abs(w[j]);
    int eu=w[i]>0?1:-1,ev=w[j]>0?1:-1;
    Poly ctab=product(remainder(w,i,j));
    Z D=ev*ctab.get(a,b);
    Z fusion_sum=0;
    for(int c=std::abs(a-b);c<=a+b;c+=2)fusion_sum+=ctab.get(c,0);
    assert(phi==2*(fusion_sum+D));

    auto [A,B]=interior_split(remainder(w,i,j));
    Poly left=product(A), base=product(B);
    Poly yp=pair_pure(base,a,b,eu,ev);
    Poly ym=pair_mixed(base,a,b,eu,ev);
    PairData z; z.phi=phi; z.D=D;
    z.P.resize(left.top+1); z.M.resize(left.top+1);
    Z budget=0,plain=0;
    for(int h=0;h<=left.top;++h) {
        z.P[h]=dot_height(left,yp,h);
        z.M[h]=dot_height(left,ym,h);
        plain+=z.P[h]+z.M[h];
        budget+=z.P[h]+(z.M[h]<0?z.M[h]:Z(0));
        ++z.prefix_count;
        if(plain<0&&z.first_plain_bad<0)z.first_plain_bad=h;
        if(budget<0&&z.first_budget_bad<0)z.first_budget_bad=h;
    }
    Z sp=std::accumulate(z.P.begin(),z.P.end(),Z(0));
    Z sm=std::accumulate(z.M.begin(),z.M.end(),Z(0));
    assert(sp+sm==phi);
    assert(sm==2*D);
    assert(budget<=phi);
    return z;
}
static std::string word_string(const std::vector<int>& w) {
    std::ostringstream s; for(int x:w)s<<x<<" "; return s.str();
}
struct Summary {
    long long words=0,pairs=0,prefixes=0,flip_words=0;
    long long budget_bad_pairs=0,plain_bad_pairs=0;
    long long budget_bad_words=0,plain_bad_words=0;
    void merge(const Summary& z) {
        words+=z.words;pairs+=z.pairs;prefixes+=z.prefixes;flip_words+=z.flip_words;
        budget_bad_pairs+=z.budget_bad_pairs;plain_bad_pairs+=z.plain_bad_pairs;
        budget_bad_words+=z.budget_bad_words;plain_bad_words+=z.plain_bad_words;
    }
}
;
struct Task { int group; Z expected; int i,j; std::vector<int> w; std::string note; };
static void receipt(const Task& job,const Z& phi) {
    auto z=evaluate_pair(job.w,job.i,job.j,phi);
    auto [A,B]=interior_split(remainder(job.w,job.i,job.j));
    std::vector<Z> available(z.P.size()); std::map<std::pair<int,int>,Z> flow;
    Z fixed_mixed=0,budget=0;
    std::cout<<"MATCHING "<<job.note<<" word="<<word_string(job.w)
             <<" pair=("<<job.w[job.i]<<","<<job.w[job.j]<<") A="<<word_string(A)
             <<" B="<<word_string(B)<<" Phi="<<phi<<"\n";
    for(int h=0;h<(int)z.P.size();++h) {
        if(z.P[h]>0)available[h]+=z.P[h];
        Z need=(z.P[h]<0?-z.P[h]:Z(0))+(z.M[h]<0?-z.M[h]:Z(0));
        for(int lower=0;lower<=h&&need>0;++lower) {
            Z take=std::min(available[lower],need);
            if(take>0){available[lower]-=take;need-=take;flow[{h,lower}]+=take;}
        }
        assert(need==0);
        if(z.M[h]>0)fixed_mixed+=z.M[h];
        budget+=z.P[h]+(z.M[h]<0?z.M[h]:Z(0));
        assert(budget>=0);
        if(z.P[h]!=0||z.M[h]!=0)
            std::cout<<"  h="<<h<<" P="<<z.P[h]<<" M="<<z.M[h]<<" Bprefix="<<budget<<"\n";
    }
    Z fixed_pure=std::accumulate(available.begin(),available.end(),Z(0));
    assert(fixed_pure+fixed_mixed==phi);
    for(auto& e:flow)if(e.second!=0)
        std::cout<<"  pair-remnants height "<<e.first.first<<" <- pure height "
                 <<e.first.second<<" count "<<e.second<<"\n";
    std::cout<<"  fixed pure="<<fixed_pure<<" fixed mixed="<<fixed_mixed
             <<" total="<<(fixed_pure+fixed_mixed)<<"\n";
}
static Poly spchar(int a,int b) {
    assert(a>=b&&b>=0); Poly p(a+b);
    for(int j=0;j<=b;++j)for(int k=0;k<=a-b;++k)
        p.at(j+k,j+a-b-k)+=1;
    return p;
}
static void genuine_module_checks() {
    long long tested=0;
    for(int a=0;a<=6;++a)for(int b=0;b<=a;++b)
    for(int c=0;c<=6;++c)for(int d=0;d<=c;++d)
    for(int r=0;r<=2;++r) {
        Poly x=spchar(a,b),y=spchar(c,d);
        for(int j=0;j<r;++j)x=signed_factor(x,-1);
        for(int j=0;j<2-r;++j)y=signed_factor(y,-1);
        Z prefix=0;
        for(int t=0;t<=std::min(x.top,y.top);++t) {
            prefix+=dot_height(x,y,t);
            assert(prefix>=0);
            ++tested;
        }
    }
    Poly v=spchar(1,1);
    for(int j=0;j<4;++j)v=signed_factor(v,-1);
    assert(v.get(0,0)==-6);
    stamp("genuine lift prefix checks="+std::to_string(tested)+
          "; exact Phi(d^4 chi_(1,1))=-6");
}

int main(int argc,char**argv) {
    int threads=24;
    for(int x=1;x<argc;++x) {
        std::string s=argv[x];
        if(s=="--threads"&&x+1<argc)threads=std::stoi(argv[++x]);
        else if(s=="-h"||s=="--help"){std::cout<<"exact interior-budget checks; --threads N\n";return 0;}
        else return 2;
    }
    omp_set_num_threads(threads);
    std::vector<Task> jobs; int group,i,j,n; Z expect;
    std::string note;
    while(std::cin>>group>>expect>>i>>j>>n) {
        std::vector<int>w(n);for(int&x:w)std::cin>>x;
        std::cin>>note;
        jobs.push_back({group,expect,i,j,w,note});
    }
    assert(jobs.size()==5751);
    stamp("begin exact cpp_int evaluator; jobs="+std::to_string(jobs.size()));
    genuine_module_checks();
    long long random_pass=0;
    for(const auto& job:jobs)if(job.group==3) {
        Poly full=product(job.w); Z phi=full.get(0,0);
        evaluate_pair(job.w,job.i,job.j,phi);
        ++random_pass;
    }
    stamp("random pair identities exact PASS count="+std::to_string(random_pass));

    Summary totals[2]; std::vector<std::string> bad_example(2);
    std::atomic<int> completed{0}; std::mutex out_lock;
    #pragma omp parallel
    {
        Summary local[2];
        #pragma omp for schedule(dynamic,1)
        for(int ix=0;ix<(int)jobs.size();++ix) {
            const auto& job=jobs[ix];
            if(job.group==3)continue;
            if(job.group==2) continue;
            int g=job.group;
            Poly full=product(job.w); Z phi=full.get(0,0);
            if(job.expected>=0)assert(phi==job.expected);
            assert(phi>=0);
            bool noflip=true,badbudget=false,badplain=false,anybudget=false;
            Summary s;s.words=1;
            for(int a=0;a<(int)job.w.size();++a)for(int b=a+1;b<(int)job.w.size();++b) {
                auto z=evaluate_pair(job.w,a,b,phi);
                ++s.pairs;s.prefixes+=z.prefix_count;
                noflip &= z.D<0;
                bool bb=z.first_budget_bad>=0, bp=z.first_plain_bad>=0;
                s.budget_bad_pairs+=bb;s.plain_bad_pairs+=bp;
                badbudget|=bb;badplain|=bp;anybudget|=!bb;
                if(bb&&bad_example[g].empty()) {
                    std::lock_guard<std::mutex> lk(out_lock);
                    if(bad_example[g].empty()) {
                        std::ostringstream o;o<<"word="<<word_string(job.w)<<" pair-index="
                            <<a<<","<<b<<" first-budget-T="<<z.first_budget_bad;
                        bad_example[g]=o.str();
                    }
                }
                if(bp&&bad_example[g].empty()) {
                    std::lock_guard<std::mutex> lk(out_lock);
                    if(bad_example[g].empty()) {
                        std::ostringstream o;o<<"word="<<word_string(job.w)<<" pair-index="
                            <<a<<","<<b<<" first-plain-T="<<z.first_plain_bad;
                        bad_example[g]=o.str();
                    }
                }
            }
            if(noflip)++s.flip_words;
            s.budget_bad_words=badbudget;
            s.plain_bad_words=badplain;
            assert(g!=0||noflip);
            local[g].merge(s);
            int d=++completed;
            if(d%400==0) {
                #pragma omp critical
                stamp("processed census/run words="+std::to_string(d));
            }
        }
        #pragma omp critical
        { totals[0].merge(local[0]); totals[1].merge(local[1]); }
    }
    std::cout<<"CENSUS words="<<totals[0].words<<" noflip="<<totals[0].flip_words
             <<" pairs="<<totals[0].pairs<<" prefixes="<<totals[0].prefixes
             <<" budget_fail_pairs="<<totals[0].budget_bad_pairs
             <<" plain_fail_pairs="<<totals[0].plain_bad_pairs
             <<" budget_fail_words="<<totals[0].budget_bad_words
             <<" plain_fail_words="<<totals[0].plain_bad_words<<"\n";
    std::cout<<"RUNS words="<<totals[1].words<<" noflip="<<totals[1].flip_words
             <<" pairs="<<totals[1].pairs<<" prefixes="<<totals[1].prefixes
             <<" budget_fail_pairs="<<totals[1].budget_bad_pairs
             <<" plain_fail_pairs="<<totals[1].plain_bad_pairs
             <<" budget_fail_words="<<totals[1].budget_bad_words
             <<" plain_fail_words="<<totals[1].plain_bad_words<<"\n";
    assert(totals[0].words==5430&&totals[0].pairs==148972);
    assert(totals[1].words==68&&totals[1].pairs==7364);
    assert(totals[0].budget_bad_pairs==0&&totals[0].plain_bad_pairs==0);
    assert(totals[1].budget_bad_pairs==0&&totals[1].plain_bad_pairs==0);
    for(int g=0;g<2;++g)if(!bad_example[g].empty())std::cout<<"FAIL "<<bad_example[g]<<"\n";

    for(const auto& job:jobs) if(job.group==2) {
        Poly f=product(job.w); Z phi=f.get(0,0);
        assert(job.expected<0||job.expected==phi);
        receipt(job,phi);
    }
    stamp("FM-STR7e independent screens PASS");
    return 0;
}
"""
# Encode each receipt label as one input token for the C++ reader.
wire = "".join(
    f"{g} {phi} {i} {j} {len(w)} " + " ".join(map(str,w)) + " " + (note or "-") + "\n"
    for g,w,phi,i,j,note in jobs
)

obj = os.memfd_create("str7e_fresh_obj",0)
exe = os.memfd_create("str7e_fresh_exe",0)
subprocess.run(
    ["g++","-O2","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],
    input=CPP,text=True,pass_fds=(obj,),check=True
)
probe = subprocess.run(
    ["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,check=True
)
link = next(shlex.split(s) for s in probe.stderr.splitlines() if "/collect2 " in s)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}","--threads",str(args.threads)],
               input=wire,text=True,pass_fds=(exe,),check=True)
