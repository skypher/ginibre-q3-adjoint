import os, subprocess, time, re, sys
from fractions import Fraction

if any(x in ("-h", "--help") for x in sys.argv[1:]):
    print("Exact in-memory FM-SEC175 scanner. Scans all admissible p and pair-free sign patterns.")
    print("Set FM_DETAIL=1 to print every no-flip pattern, TopPair, and exact ratio.")
    raise SystemExit(0)

CPP = r'''#include <gmpxx.h>
#include <bits/stdc++.h>
using namespace std;
using Z = mpz_class;

Z coeff(const vector<vector<Z>>& r, const vector<int>& w,
        const vector<int>& pa, int S, int k) {
    if (k < 0 || k > w[S] || ((k - pa[S]) & 1)) return 0;
    int q = (k - pa[S]) / 2;
    return q < 0 || q >= (int)r[S].size() ? Z(0) : r[S][q];
}
void make_rows(const vector<int>& ns, vector<vector<Z>>& r,
               vector<int>& w, vector<int>& pa) {
    int n = ns.size(), N = 1 << n;
    r.assign(N, {}); w.assign(N, 0); pa.assign(N, 0);
    r[0] = {Z(1)};
    for (int S = 1; S < N; ++S) {
        int b = __builtin_ctz((unsigned)S), old = S ^ (1 << b), a = ns[b];
        w[S] = w[old] + a; pa[S] = w[S] & 1;
        vector<Z> pref(r[old].size() + 1);
        for (size_t j = 0; j < r[old].size(); ++j)
            pref[j + 1] = pref[j] + r[old][j];
        int len = (w[S] - pa[S]) / 2 + 1;
        r[S].resize(len);
        for (int j = 0; j < len; ++j) {
            int k = pa[S] + 2*j, lo = abs(k-a), hi = min(w[old], k+a);
            int qlo = lo <= pa[old] ? 0 : (lo-pa[old]+1)/2;
            int qhi = (hi-pa[old])/2;
            if (qhi >= 0 && qlo < (int)r[old].size()) {
                qhi = min(qhi, (int)r[old].size()-1);
                if (qlo <= qhi) r[S][j] = pref[qhi+1] - pref[qlo];
            }
        }
    }
}
void fwht(vector<Z>& v) {
    for (int h=1; h<(int)v.size(); h*=2)
        for (int a=0; a<(int)v.size(); a+=2*h)
            for (int j=0; j<h; ++j) {
                Z x=v[a+j], y=v[a+j+h];
                v[a+j]=x+y; v[a+j+h]=x-y;
            }
}
vector<Z> g_values(int full, int p, const vector<vector<Z>>& r,
                   const vector<int>& w, const vector<int>& pa) {
    vector<Z> v(full+1);
    for (int S=0; S<=full; ++S)
        v[S] = coeff(r,w,pa,S,0) * coeff(r,w,pa,full^S,p);
    fwht(v);
    return v;
}
string zs(const Z& x) { return x.get_str(); }

int main(int argc, char** argv) {
    if (argc>1 && (string(argv[1])=="-h" || string(argv[1])=="--help")) {
        cout << "Usage: scanner n1:m1 n2:m2 ... (GMP exact)." << endl;
        return 0;
    }
    vector<int> cl, mult, ns;
    for (int i=1; i<argc; ++i) {
        int a,b;
        if (sscanf(argv[i], "%d:%d", &a,&b)!=2 || a<1 || b<1) return 2;
        cl.push_back(a); mult.push_back(b);
    }
    int C=cl.size(), k=accumulate(mult.begin(),mult.end(),0);
    int W=0, mx=0, c3=0;
    for (int i=0; i<C; ++i) {
        W += cl[i]*mult[i]; mx=max(mx,cl[i]);
        if (cl[i]>=3) c3 += mult[i];
        for (int j=0; j<mult[i]; ++j) ns.push_back(cl[i]);
    }
    if (k<2 || k>16) return 2;

    vector<vector<Z>> rows;
    vector<int> wt, parity;
    make_rows(ns, rows, wt, parity);
    int full=(1<<k)-1;
    vector<int> classbits(C);
    int bit=0;
    for (int i=0; i<C; ++i)
        for (int j=0; j<mult[i]; ++j) classbits[i] |= 1<<bit++;

    int ti=-1,tj=-1,best=-1,bestmax=-1;
    for (int i=0; i<k; ++i) for (int j=i+1; j<k; ++j)
        if ((ns[i]+ns[j])%2==0) {
            int sum=ns[i]+ns[j], mxpair=max(ns[i],ns[j]);
            if (sum>best || (sum==best && mxpair>bestmax)) {
                best=sum; bestmax=mxpair; ti=i; tj=j;
            }
        }
    if (ti<0) return 3;

    vector<int> keep, child_labels;
    for (int i=0; i<k; ++i) if (i!=ti && i!=tj) {
        keep.push_back(i); child_labels.push_back(ns[i]);
    }
    vector<vector<Z>> child_rows;
    vector<int> child_wt, child_parity;
    make_rows(child_labels,child_rows,child_wt,child_parity);
    int child_full=(1<<(k-2))-1;

    int lo=max(6,mx);
    if ((lo&1)!=(W&1)) ++lo;
    int hi=W-2*max(8,mx);
    if ((hi&1)!=(W&1)) --hi;

    long long p_runs=0, patterns=0, noflip=0, tp_fail=0, ft_fail=0, nonpositive=0;
    for (int p=lo; p<=hi; p+=2) {
        if ((W-p)/2<8 || mx>(W-p)/2 || c3<2) continue;
        ++p_runs;
        vector<Z> g=g_values(full,p,rows,wt,parity);
        vector<Z> gc=g_values(child_full,p,child_rows,child_wt,child_parity);
        long long pn=0,nf=0,fl=0,ft=0,np=0;
        Z max_num=-1,max_den=1;
        string max_record="none";

        for (int s=0; s<(1<<C); ++s) {
            int neg=0; unsigned mask=0;
            for (int i=0; i<C; ++i) if (s>>i&1) {
                neg+=mult[i]; mask|=classbits[i];
            }
            int sigma=(neg&1)?-1:1;
            bool pair=false;
            for (int i=0; i<C; ++i)
                if (cl[i]==p && ((s>>i&1)?-1:1)!=sigma) pair=true;
            if (pair) continue;
            ++pn;

            bool descent=false;
            for (int i=0; i<k && !descent; ++i)
                for (int j=i+1; j<k; ++j)
                    if (g[mask]-g[mask^(1u<<i)^(1u<<j)]>=0) {
                        descent=true; break;
                    }
            if (!descent)
                for (int i=0; i<k; ++i)
                    if (g[mask]-g[mask^(1u<<i)]>=0) {
                        descent=true; break;
                    }
            if (descent) continue;
            ++nf;

            unsigned cm=0; int z=0;
            for (int j=0; j<k; ++j) if (j!=ti && j!=tj) {
                if (mask>>j&1) cm|=1u<<z;
                ++z;
            }
            Z parent=g[mask], child=gc[cm];
            if (child>parent) ++fl;
            if (16*child>parent) ++ft;

            string minus;
            for (int i=0; i<C; ++i) if (s>>i&1) {
                if (!minus.empty()) minus+=",";
                minus+=to_string(cl[i])+"^"+to_string(mult[i]);
            }
            int ta=ns[ti], tb=ns[tj];
            if (mask>>ti&1) ta=-ta;
            if (mask>>tj&1) tb=-tb;
            if (parent<=0) {
                ++np;
                if (getenv("FM_DETAIL"))
                    cout << "NOFLIP p="<<p<<" minus="<<(minus.empty()?"none":minus)
                         <<" sigma_p="<<sigma*p<<" TopPair=("<<ta<<","<<tb
                         <<") parent="<<zs(parent)<<" child="<<zs(child)
                         <<" ratio=undefined"<<endl;
                continue;
            }
            if (child*max_den>max_num*parent) {
                max_num=child; max_den=parent;
                max_record="minus="+(minus.empty()?string("none"):minus)
                    +" sigma_p="+to_string(sigma*p)
                    +" TopPair=("+to_string(ta)+","+to_string(tb)+")"
                    +" parent="+zs(parent)+" child="+zs(child);
            }
            if (getenv("FM_DETAIL"))
                cout << "NOFLIP p="<<p<<" minus="<<(minus.empty()?"none":minus)
                     <<" sigma_p="<<sigma*p<<" TopPair=("<<ta<<","<<tb
                     <<") parent="<<zs(parent)<<" child="<<zs(child)
                     <<" ratio="<<zs(child)<<"/"<<zs(parent)<<endl;
        }
        patterns+=pn; noflip+=nf; tp_fail+=fl; ft_fail+=ft; nonpositive+=np;
        if (nf)
            cout << "P_RESULT p="<<p<<" patterns="<<pn<<" noflip="<<nf
                 <<" TPfail="<<fl<<" FTprimefail="<<ft
                 <<" nonpositive_parent="<<np<<" maxratio="
                 <<(max_num<0?string("none"):zs(max_num)+"/"+zs(max_den))
                 <<" "<<max_record<<endl;
    }
    cout << "SUMMARY p_runs="<<p_runs<<" patterns="<<patterns
         <<" noflip="<<noflip<<" TPfail="<<tp_fail
         <<" FTprimefail="<<ft_fail<<" nonpositive_parent="<<nonpositive<<endl;
}
'''

fd=os.memfd_create("fm_sec175_scan")
os.set_inheritable(fd,True)
env=dict(os.environ,TMPDIR="/dev/shm",OMP_NUM_THREADS="4")
build=subprocess.run(
    ["g++","-pipe","-std=c++17","-O3","-x","c++","-o",f"/proc/self/fd/{fd}","-","-lgmpxx","-lgmp"],
    input=CPP.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,
    pass_fds=(fd,),env=env)
if build.returncode: raise SystemExit(build.stderr.decode())
os.fchmod(fd,0o700)
smoke=subprocess.run([f"/proc/self/fd/{fd}","1:2","2:1","3:4","4:3","5:4"],
    pass_fds=(fd,),env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
assert smoke.returncode==0 and "parent=4075371 child=160741" in smoke.stdout,smoke.stdout
print("compiled; smoke check: exact witness 4075371 -> 160741",flush=True)

def band(W): return "60-99" if W<100 else ("100-149" if W<150 else "150-200")
runs={}
for k in range(8,13):
    for step in (1,2):
        for s in range(1,201):
            labels=tuple(range(s,s+step*k,step)); W=sum(labels); mx=labels[-1]
            lo=max(6,mx); lo+=((W-lo)&1)
            hi=W-2*max(8,mx); hi-=((W-hi)&1)
            if 60<=W<=200 and hi>=lo and sum(x>=3 for x in labels)>=2:
                runs[tuple((x,1) for x in labels)]=(k,step,s,W)

mixed={}
for ones in (1,2):
    for nr in range(4,10):
        if not 8<=nr+ones+2<=12: continue
        for step in (1,2):
            for s in range(4,201):
                d={1:ones,2:1,3:1}
                for x in range(s,s+step*nr,step): d[x]=d.get(x,0)+1
                W=sum(x*y for x,y in d.items()); mx=max(d)
                lo=max(6,mx); lo+=((W-lo)&1)
                hi=W-2*max(8,mx); hi-=((W-hi)&1)
                if 60<=W<=200 and hi>=lo and sum(y for x,y in d.items() if x>=3)>=2:
                    mixed[tuple(sorted(d.items()))]=(ones,nr,step,s,W)

print("profiles",len(runs),len(mixed),flush=True)
summary={name:{b:[0,0,0,0,0,Fraction(-1),None]
               for b in ("60-99","100-149","150-200")}
         for name in ("runs","mixed")}
for name, profiles in (("runs",runs),("mixed",mixed)):
    start=time.monotonic()
    for i,(items,tag) in enumerate(sorted(profiles.items(),key=lambda q:(q[1][-1],q[0])),1):
        if i==1 or i%10==0 or i==len(profiles):
            print("PROGRESS",time.strftime("%H:%M:%S"),name,i,"/",len(profiles),tag,flush=True)
        z=subprocess.run([f"/proc/self/fd/{fd}",*[f"{x}:{y}" for x,y in items]],
            pass_fds=(fd,),env=env,text=True,stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,timeout=180)
        if z.returncode: raise RuntimeError(z.stderr+"\n"+z.stdout)
        q=summary[name][band(tag[-1])]
        for line in z.stdout.splitlines():
            if line.startswith("SUMMARY "):
                mm=re.search(r"p_runs=(\d+) patterns=(\d+) noflip=(\d+) TPfail=(\d+) FTprimefail=(\d+) nonpositive_parent=(\d+)",line)
                assert mm,line
                vals=tuple(map(int,mm.groups()))
                for j,v in enumerate(vals[:5]): q[j]+=v
            elif line.startswith("P_RESULT "):
                mm=re.search(r"p=(\d+).*maxratio=(\S+) (.*)",line); assert mm,line
                p,rat,rec=mm.groups()
                if rat!="none":
                    a,b=rat.split("/"); rat=Fraction(int(a),int(b))
                    if rat>q[5]: q[5]=rat; q[6]=(tag,int(p),rec)
            elif line.startswith("NOFLIP "):
                print(name,line,flush=True)
    print("FAMILY_DONE",name,"seconds",round(time.monotonic()-start,2),flush=True)

for name in ("runs","mixed"):
    print("BANDS",name)
    for b,q in summary[name].items():
        print(b,{"p_runs":q[0],"patterns":q[1],"noflip":q[2],
                 "TopPair_fail":q[3],"FTprime_fail":q[4],
                 "maxratio":str(q[5]),"record":q[6]})
    print("TOTAL",name,{key:sum(q[j] for q in summary[name].values())
          for key,j in (("p_runs",0),("patterns",1),("noflip",2),
                        ("TopPair_fail",3),("FTprime_fail",4))})
PY
