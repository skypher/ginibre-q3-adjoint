import ast, os, subprocess, sys
from pathlib import Path

base=Path('ginibre_q3/character_ring_iter/fm39/sec160_S_exhaustive.py').read_text()
tree=ast.parse(base)
src=ast.literal_eval(next(n.value for n in tree.body if isinstance(n,ast.Assign)
    and any(isinstance(t,ast.Name) and t.id=='src' for t in n.targets)))

src=src.replace('''    BranchStats branch[3];
    Failure failure;
};''','''    BranchStats branch[3];
    Failure failure;
    long long residual_rows=0, d_fail=0;
    array<long long,64> d_hist{};
    array<long long,4096> mask_hist{};
    array<long long,12> cand_hits{};
};''')

helper=r'''
vector<Rem> candidateSet(const Rec& r, const vector<Rem>& v) {
    vector<Rem> z(12, Rem{-1,-1});
    auto put=[&](int i, Rem R) {
        for (auto q : v) if (q.n==R.n && q.m==R.m) { z[i]=R; return; }
    };
    put(0, ruleW(r,v).R); put(1, ruleM(r,v).R);
    int lo=0, hi=0, bn=0, bw=-1, mx=0;
    for(int n=1;n<=L;++n) if(r.c[n]) {
        mx=n;
        int f=abs((int)r.c[n]), w=f*n;
        if(f>=2 && (!lo || n<lo)) lo=n;
        if(f>=2 && n>hi) hi=n;
        if(w>bw || (w==bw && n>bn)) bw=w,bn=n;
    }
    if(lo) put(2,{lo,lo});
    if(hi) put(3,{hi,hi});
    if(bn && abs((int)r.c[bn])>=2) put(4,{bn,bn});
    int mf=-1,mn=-1;
    for(int n=1;n<=L;++n) if(r.c[n]) {
        int f=abs((int)r.c[n]);
        if(f>mf || (f==mf && n>mn)) mf=f,mn=n;
    }
    if(mn>=0 && mf>=2) put(5,{mn,mn});
    Rem pl{-1,-1},pg{-1,-1};
    for(auto q:v) if(q.m) {
        if(pl.n<0 || q.n<pl.n || (q.n==pl.n && q.m<pl.m)) pl=q;
        if(pg.n<0 || q.n>pg.n || (q.n==pg.n && q.m>pg.m)) pg=q;
    }
    if(pl.n>=0) put(6,pl);
    if(pg.n>=0) put(7,pg);
    if(mx && abs((int)r.c[mx])==1)
        for(int n=mx-1;n>=1;--n)
            if(r.c[n] && ((n^mx)&1)==0) { put(8,{n,mx}); break; }
    Rem es{-1,-1},el{-1,-1};
    for(auto q:v) if(!q.m) {
        if(es.n<0 || q.n<es.n) es=q;
        if(q.n>el.n) el=q;
    }
    if(es.n>=0) put(9,es);
    if(el.n>=0) put(10,el);
    if(mx) {
        if(abs((int)r.c[mx])>=2) put(11,{mx,mx});
        else for(int n=mx-1;n>=1;--n)
            if(r.c[n] && ((n^mx)&1)==0) { put(11,{n,mx}); break; }
    }
    return z;
}
'''
src=src.replace('string show(const Rec& r) {',helper+'\nstring show(const Rec& r) {')

src=src.replace('''                auto v = allowed(r);
                if (v.empty()) continue;''','''                auto v = allowed(r);
                if (v.empty()) continue;
                auto csel = candidateSet(r,v);''')

src=src.replace('''                    I parent = gp(r, p), childM = gp(childRec, p), childW = gp(childRec2, p);
                    I child = childM + childW - parent;''','''                    I parent = gp(r, p), childM = gp(childRec, p), childW = gp(childRec2, p);
                    I child = childM + childW - parent;
                    int dd0=(r.W-p)/2, kc0=0;
                    for(int nn=3;nn<=L;++nn) kc0+=abs((int)r.c[nn]);
                    if(p>=6 && dd0>=8 && r.mx<=dd0 && kc0>=2) {
                        ++s.residual_rows; int nm=0;
                        for(auto rr:v) {
                            auto ci=lookup.find(childKey(r,rr));
                            if(ci!=lookup.end() && gp(recs[ci->second],p)<=parent) ++nm;
                        }
                        if(!nm) ++s.d_fail;
                        ++s.d_hist[min(nm,63)];
                        int mask=0;
                        for(int jj=0;jj<12;++jj) if(csel[jj].n>=0) {
                            auto ci=lookup.find(childKey(r,csel[jj]));
                            if(ci!=lookup.end() && gp(recs[ci->second],p)<=parent) mask|=(1<<jj);
                        }
                        ++s.mask_hist[mask];
                        for(int jj=0;jj<12;++jj) if(mask&(1<<jj)) ++s.cand_hits[jj];
                    }''')

src=src.replace('''    Failure firstUnproved;
    BranchStats total[3];''','''    Failure firstUnproved;
    long long residual_rows=0, d_fail=0;
    array<long long,64> d_hist{};
    array<long long,4096> mask_hist{};
    array<long long,12> cand_hits{};
    BranchStats total[3];''')

src=src.replace('''                eligible += s.eligible;
                tests += s.tests;''','''                eligible += s.eligible;
                tests += s.tests;
                residual_rows+=s.residual_rows; d_fail+=s.d_fail;
                for(int q=0;q<64;++q) d_hist[q]+=s.d_hist[q];
                for(int q=0;q<4096;++q) mask_hist[q]+=s.mask_hist[q];
                for(int q=0;q<12;++q) cand_hits[q]+=s.cand_hits[q];''')

old='''    cout << flush;
}
'''
new=r'''    const char* cn[12]={"RuleW","RuleM","dup_small","dup_large","dup_heaviest","dup_most_frequent","pair_smallest","pair_largest","max_singleton_largest_same_parity","even_small","even_large","top_two_compatible"};
    cout << "D residual_rows=" << residual_rows << " D_failures=" << d_fail << "\n";
    for(int j=1;j<64;++j) if(d_hist[j]) cout << "D_monotone_count " << j << " rows " << d_hist[j] << "\n";
    for(int j=0;j<12;++j) cout << "candidate " << j << " " << cn[j] << " success_rows=" << cand_hits[j] << "\n";
    int best=13,nsol=0; vector<int> sol;
    for(int m=1;m<4096;++m) {
        bool hit=true;
        for(int k=0;k<4096;++k) if(mask_hist[k] && !(m&k)) {hit=false;break;}
        if(hit) {
            int n=__builtin_popcount((unsigned)m);
            if(n<best) best=n,nsol=0,sol.clear();
            if(n==best) {++nsol;sol.push_back(m);}
        }
    }
    cout << "minimum candidate set size=" << best << " count=" << nsol << "\n";
    for(int m:sol) {
        cout << "SET";
        for(int j=0;j<12;++j) if(m&(1<<j)) cout << " " << j << ":" << cn[j];
        cout << "\n";
    }
    int distinct=0;
    for(int k=0;k<4096;++k) if(mask_hist[k]) ++distinct;
    cout << "distinct_success_masks=" << distinct << "\n";
    for(int k=0;k<4096;++k) if(mask_hist[k]) cout << "mask " << k << " count " << mask_hist[k] << "\n";
    cout << flush;
}
'''
if old not in src: raise SystemExit('missing final marker')
src=src.replace(old,new,1)

fd=os.memfd_create('fm_sec163_exact',0)
os.set_inheritable(fd,True)
env=dict(os.environ)
env['TMPDIR']='/dev/shm'
env['OMP_NUM_THREADS']='32'
build=subprocess.run(
    ['g++','-pipe','-std=c++17','-O3','-fopenmp','-x','c++','-o',
     f'/proc/self/fd/{fd}','-'],
    input=src.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,
    pass_fds=(fd,),env=env)
if build.returncode:
    print(build.stderr.decode(),file=sys.stderr)
    raise SystemExit(build.returncode)
os.fchmod(fd,0o700)
print('compiled exact FM-SEC163 sweep; MAXW=40; OpenMP=32',flush=True)
run=subprocess.run([f'/proc/self/fd/{fd}'],pass_fds=(fd,),env=env)
raise SystemExit(run.returncode)
