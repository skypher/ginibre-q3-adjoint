import ast
import os
import pathlib
import subprocess
import sys

path = pathlib.Path("ginibre_q3/character_ring_iter/fm39/sec160_WM_exhaustive48.py")
tree = ast.parse(path.read_text())
src = None
for node in tree.body:
    if isinstance(node, ast.Assign) and any(
        isinstance(t, ast.Name) and t.id == "src" for t in node.targets
    ):
        src = ast.literal_eval(node.value)
        break
assert src is not None
src = src.replace(
    "static const int MAXW = 48, L = 48;",
    "static const int MAXW = 44, L = 44;"
)
pos = src.index("int main() {")
src = src[:pos] + r"""
struct SweepCounts {
    long long residual_backgrounds=0, residual_pairs=0, comparisons=0;
    long long failures=0, no_removal=0, key_misses=0, parent_nonpositive=0;
    long long argmax_terms=0, tied_cases=0;
    long long type[9]{};
    long long sel_defined[5]{}, sel_exact[5]{}, sel_half[5]{};
};
void addCounts(SweepCounts& a, const SweepCounts& b) {
    a.residual_backgrounds+=b.residual_backgrounds;
    a.residual_pairs+=b.residual_pairs;
    a.comparisons+=b.comparisons;
    a.failures+=b.failures;
    a.no_removal+=b.no_removal;
    a.key_misses+=b.key_misses;
    a.parent_nonpositive+=b.parent_nonpositive;
    a.argmax_terms+=b.argmax_terms;
    a.tied_cases+=b.tied_cases;
    for(int i=0;i<9;++i) a.type[i]+=b.type[i];
    for(int i=0;i<5;++i) {
        a.sel_defined[i]+=b.sel_defined[i];
        a.sel_exact[i]+=b.sel_exact[i];
        a.sel_half[i]+=b.sel_half[i];
    }
}
int category(const Rec& r, Rem x) {
    if(!x.m) return 0;
    if(x.n==x.m) {
        int smallestDup=0, heavy=0;
        long long bestWeight=-1;
        for(int n=1;n<=L;++n)
            if(abs((int)r.c[n])>=2 && !smallestDup) smallestDup=n;
        for(int n=1;n<=L;++n) if(r.c[n]) {
            long long wt=(long long)abs((int)r.c[n])*n;
            if(wt>bestWeight || (wt==bestWeight && n>heavy)) {
                bestWeight=wt; heavy=n;
            }
        }
        if(x.n==smallestDup) return 1;
        if(x.n==heavy) return 2;
        if(x.n==r.mx) return 3;
        return 4;
    }
    int first=0, second=0;
    for(int n=1;n<=L;++n) if(r.c[n]) {
        if(!first) {
            first=n;
            if(abs((int)r.c[n])>=2) { second=n; break; }
        } else { second=n; break; }
    }
    if(first && second && x.n==first && x.m==second) return 5;
    int maxSame=0;
    for(int n=1;n<r.mx;++n)
        if(r.c[n] && (n&1)==(r.mx&1)) maxSame=n;
    if(maxSame && x.n==maxSame && x.m==r.mx) return 6;
    int heavy=0;
    long long bestWeight=-1;
    for(int n=1;n<=L;++n) if(r.c[n]) {
        long long wt=(long long)abs((int)r.c[n])*n;
        if(wt>bestWeight || (wt==bestWeight && n>heavy)) {
            bestWeight=wt; heavy=n;
        }
    }
    if(heavy && abs((int)r.c[heavy])==1) {
        int partner=0;
        for(int n=1;n<=L;++n)
            if(n!=heavy && r.c[n] && (n&1)==(heavy&1)) partner=n;
        if(partner && x.n==min(heavy,partner) && x.m==max(heavy,partner))
            return 7;
    }
    return 8;
}
int remIndex(const vector<Rem>& v, Rem x) {
    for(int i=0;i<(int)v.size();++i)
        if(v[i].n==x.n && v[i].m==x.m) return i;
    return -1;
}
int main() {
    U radix=1;
    for(int n=1;n<=L;++n) {
        mult[n]=radix;
        radix*=U(2*(MAXW/n)+1);
    }
    if(radix >= (U(1)<<127)) return 2;
    recs.reserve(15200000);
    auto start=chrono::steady_clock::now();
    generate(1,0,0,0,0);
    sort(recs.begin(),recs.end(),[](const Rec&a,const Rec&b) {
        return a.W!=b.W ? a.W<b.W : a.key<b.key;
    });
    cout<<"generated_profiles="<<recs.size()<<" max_weight="<<MAXW<<"\n"<<flush;
    unordered_map<U,size_t,HashU> index;
    index.reserve((size_t)(recs.size()*1.3));
    for(size_t i=0;i<recs.size();++i) index[recs[i].key]=i;
    const auto& lookup=index;
    SweepCounts total;
    const char* snames[5]={
        "Rule_W","Rule_M","Rule_1prime",
        "minimum_total_removal_weight","minimum_weight_same_parity_pair"
    };
    int nt=max(1,omp_get_max_threads());
    const size_t CHUNK=25000;
    for(int W=0;W<=MAXW;++W) {
        auto lower=[&](int w) {
            return lower_bound(recs.begin(),recs.end(),w,
                [](const Rec&r,int x){return r.W<x;})-recs.begin();
        };
        size_t lo=lower(W), hi=lower(W+1);
        for(size_t base=lo;base<hi;base+=CHUNK) {
            size_t end=min(hi,base+CHUNK);
            vector<SweepCounts> local(nt);
            #pragma omp parallel for num_threads(nt) schedule(dynamic,8)
            for(long long ii=(long long)base;ii<(long long)end;++ii) {
                Rec& r=recs[ii];
                r.g=calculate(r);
                int coreCount=0;
                for(int n=3;n<=L;++n) coreCount+=abs((int)r.c[n]);
                if(coreCount<2) continue;
                int pmin=max(6,r.mx), pmax=min(r.W-16,r.W-2*r.mx);
                if((pmin&1)!=(r.W&1)) ++pmin;
                if(pmin>pmax) continue;
                SweepCounts& s=local[omp_get_thread_num()];
                ++s.residual_backgrounds;
                vector<Rem> v=allowed(r);
                if(v.empty()) {
                    s.no_removal+=((pmax-pmin)/2+1);
                    continue;
                }
                vector<const Rec*> childRec(v.size(),nullptr);
                bool missing=false;
                for(int k=0;k<(int)v.size();++k) {
                    auto it=lookup.find(childKey(r,v[k]));
                    if(it==lookup.end()) { missing=true; break; }
                    childRec[k]=&recs[it->second];
                }
                if(missing) { ++s.key_misses; continue; }

                Pick rw=ruleW(r,v), rm=ruleM(r,v), r1=rule1prime(r,v);
                Rem minTotal=v[0], minPair=v[0];
                bool havePair=false;
                auto remWeight=[](Rem x){return x.m?x.n+x.m:x.n;};
                for(Rem x:v) {
                    if(remWeight(x)<remWeight(minTotal) ||
                       (remWeight(x)==remWeight(minTotal) &&
                        (x.n<minTotal.n ||
                         (x.n==minTotal.n && x.m<minTotal.m))))
                        minTotal=x;
                    if(x.m && (!havePair || x.n+x.m<minPair.n+minPair.m ||
                       (x.n+x.m==minPair.n+minPair.m &&
                        (x.n<minPair.n || (x.n==minPair.n && x.m<minPair.m))))) {
                        minPair=x; havePair=true;
                    }
                }
                if(!havePair) minPair=smallestEven(v);
                int si[5]={
                    remIndex(v,rw.R), remIndex(v,rm.R), remIndex(v,r1.R),
                    remIndex(v,minTotal), remIndex(v,minPair)
                };
                for(int p=pmin;p<=pmax;p+=2) {
                    ++s.residual_pairs;
                    I parent=gp(r,p);
                    if(parent<=0) ++s.parent_nonpositive;
                    vector<I> drops(v.size());
                    I best=0; bool first=true; int ties=0;
                    for(int k=0;k<(int)v.size();++k) {
                        I d=parent-gp(*childRec[k],p);
                        drops[k]=d;
                        ++s.comparisons;
                        if(first || d>best) { best=d; ties=1; first=false; }
                        else if(d==best) ++ties;
                    }
                    if(best<0) ++s.failures;
                    s.argmax_terms+=ties;
                    if(ties>1) ++s.tied_cases;
                    for(int k=0;k<(int)v.size();++k)
                        if(drops[k]==best) ++s.type[category(r,v[k])];
                    for(int q=0;q<5;++q) if(si[q]>=0) {
                        ++s.sel_defined[q];
                        if(drops[si[q]]==best) ++s.sel_exact[q];
                        if(2*drops[si[q]]>=best) ++s.sel_half[q];
                    }
                }
            }
            for(const auto& x:local) addCounts(total,x);
            double sec=chrono::duration<double>(
                chrono::steady_clock::now()-start).count();
            cout<<"progress W="<<W<<" profiles="<<(end-lo)<<"/"<<(hi-lo)
                <<" residual_pairs="<<total.residual_pairs
                <<" removal_checks="<<total.comparisons
                <<" argmax_terms="<<total.argmax_terms
                <<" elapsed_s="<<fixed<<setprecision(1)<<sec<<"\n"<<flush;
        }
        if(W==MAXW || W==40) {
            cout<<"SUMMARY upto_W="<<W
                <<" pair_free_profiles="<<lower(W+1)
                <<" residual_backgrounds="<<total.residual_backgrounds
                <<" residual_pairs="<<total.residual_pairs
                <<" removal_comparisons="<<total.comparisons
                <<" D_failures="<<total.failures
                <<" no_removal="<<total.no_removal
                <<" key_misses="<<total.key_misses
                <<" parent_nonpositive="<<total.parent_nonpositive
                <<" argmax_terms="<<total.argmax_terms
                <<" tied_cases="<<total.tied_cases<<"\n";
            const char* types[9]={
                "single_even","duplicate_smallest_class",
                "duplicate_heaviest_class","duplicate_max_class",
                "duplicate_other_class","pair_two_smallest_factors",
                "pair_max_with_next_same_parity",
                "pair_heaviest_singleton_with_partner","pair_other"
            };
            for(int i=0;i<9;++i)
                cout<<"argmax_type "<<types[i]<<"="<<total.type[i]<<"\n";
            for(int i=0;i<5;++i)
                cout<<"selector "<<snames[i]
                    <<" defined="<<total.sel_defined[i]
                    <<" exact_argmax="<<total.sel_exact[i]
                    <<" within_half="<<total.sel_half[i]<<"\n";
            cout<<"no_residual_D_failure_through="<<W<<"\n"<<flush;
        }
    }
}
"""
fd = os.memfd_create("sec165")
os.set_inheritable(fd, True)
env = os.environ.copy()
env["TMPDIR"] = "/dev/shm"
env.setdefault("OMP_NUM_THREADS", "8")
build = subprocess.run(
    ["g++","-pipe","-std=c++17","-O3","-fopenmp","-x","c++",
     "-o",f"/proc/self/fd/{fd}","-"],
    input=src.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    pass_fds=(fd,), env=env
)
if build.returncode:
    sys.stderr.write(build.stderr.decode())
    raise SystemExit(build.returncode)
os.fchmod(fd,0o755)
run=subprocess.run([f"/proc/self/fd/{fd}"],pass_fds=(fd,),env=env)
raise SystemExit(run.returncode)
