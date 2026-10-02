import os, re, subprocess, time
from fractions import Fraction

CPP = r'''
#include <gmpxx.h>
#include <bits/stdc++.h>
using namespace std; using Z=mpz_class;

static vector<Z> calcphi(const vector<int>& ns) {
  int n=(int)ns.size(), N=1<<n, full=N-1;
  vector<int> wt(N), parity(N);
  vector<vector<Z>> rows(N); rows[0]={Z(1)};
  for(int S=1;S<N;++S) {
    int bit=__builtin_ctz((unsigned)S), old=S^(1<<bit), a=ns[bit];
    wt[S]=wt[old]+a; parity[S]=wt[S]&1;
    const vector<Z>& r=rows[old]; int po=parity[old], oldmax=wt[old];
    vector<Z> pref(r.size()+1);
    for(size_t q=0;q<r.size();++q) pref[q+1]=pref[q]+r[q];
    int len=(wt[S]-parity[S])/2+1; rows[S].resize(len);
    for(int q=0;q<len;++q) {
      int k=parity[S]+2*q, lo=abs(k-a), hi=min(oldmax,k+a);
      int qlo=(lo<=po)?0:(lo-po+1)/2, qhi=(hi-po)/2;
      if(qhi>=0 && qlo<(int)r.size()) {
        qhi=min(qhi,(int)r.size()-1);
        if(qlo<=qhi) rows[S][q]=pref[qhi+1]-pref[qlo];
      }
    }
  }
  vector<Z> f(N);
  for(int S=0;S<N;++S) {
    Z a=parity[S]?Z(0):rows[S][0];
    Z b=parity[full^S]?Z(0):rows[full^S][0];
    f[S]=a*b;
  }
  vector<vector<Z>>().swap(rows);
  for(int h=1;h<N;h*=2)
    for(int a=0;a<N;a+=2*h)
      for(int j=0;j<h;++j) {
        Z x=f[a+j], y=f[a+j+h];
        f[a+j]=x+y; f[a+j+h]=x-y;
      }
  return f;
}

static pair<int,int> topPair(const vector<int>& B) {
  int u=-1,v=-1,best=-1,bestmax=-1;
  for(int i=0;i<(int)B.size();++i)
    for(int j=i+1;j<(int)B.size();++j)
      if((B[i]+B[j])%2==0) {
        int sc=B[i]+B[j], mx=max(B[i],B[j]);
        if(sc>best || (sc==best && mx>bestmax)) {
          best=sc; bestmax=mx; u=i; v=j;
        }
      }
  return {u,v};
}

int main(int argc,char**argv) {
  if(argc<3) return 2;
  int p=atoi(argv[1]), W=0, mx=0, c3=0;
  vector<int> labels,mults,B;
  for(int z=2;z<argc;++z) {
    int n,m; sscanf(argv[z],"%d:%d",&n,&m);
    labels.push_back(n); mults.push_back(m); W+=n*m; mx=max(mx,n);
    if(n>=3)c3+=m;
    for(int t=0;t<m;++t) B.push_back(n);
  }
  int delta=(W-p)/2;
  bool residual=p>=6 && p>=mx && (W-p)%2==0 &&
                delta>=8 && mx<=delta && c3>=2;
  printf("profile W=%d p=%d delta=%d classes=%zu factors=%zu residual=%s\n",
         W,p,delta,labels.size(),B.size()+1,residual?"yes":"NO");
  if(!residual)return 0;

  vector<int> all=B; all.push_back(p);
  int nb=(int)B.size(), n=(int)all.size();
  vector<Z> f=calcphi(all);
  auto [tu,tv]=topPair(B);
  long long patterns=0,noflip=0,tpfail=0;
  Z maxnum=0,maxden=1,maxparent=0,maxchild=0;
  int maxA=0,maxB=0;

  for(int C=0;C<(1<<(int)labels.size());++C) {
    int neg=0;
    for(int cl=0;cl<(int)labels.size();++cl)
      if(C>>cl&1) neg+=mults[cl];
    if(neg<2 || neg>6)continue;

    unsigned bm=0; int bitpos=0;
    for(int cl=0;cl<(int)labels.size();++cl)
      for(int t=0;t<mults[cl];++t,++bitpos)
        if(C>>cl&1) bm|=1u<<bitpos;

    int sp=(neg&1)?-1:1; bool pair=false;
    for(int i=0;i<nb;++i)
      if(B[i]==p && (((bm>>i)&1)?-1:1)!=sp) pair=true;
    if(pair)continue;

    ++patterns;
    unsigned mask=bm; if(sp<0)mask|=1u<<nb;
    bool descent=false;
    for(int i=0;i<n && !descent;++i)
      for(int j=i+1;j<n;++j)
        if(f[mask]-f[mask^(1u<<i)^(1u<<j)]>=0) {
          descent=true; break;
        }
    if(descent)continue;
    ++noflip;
    if(tu<0) {++tpfail;continue;}

    int negRemoved=((bm>>tu)&1)+((bm>>tv)&1);
    int negChild=neg-negRemoved;
    vector<int> child; unsigned cm=0; int ci=0;
    for(int i=0;i<nb;++i) if(i!=tu && i!=tv) {
      child.push_back(B[i]);
      if((bm>>i)&1)cm|=1u<<ci;
      ++ci;
    }
    child.push_back(p); if(negChild&1)cm|=1u<<ci;

    vector<Z> fc=calcphi(child);
    Z pg=f[mask]/2, cg=fc[cm]/2;
    if(cg>pg)++tpfail;
    if(cg*maxden>maxnum*pg) {
      maxnum=cg; maxden=pg; maxparent=pg; maxchild=cg;
      maxA=((bm>>tu)&1)?-B[tu]:B[tu];
      maxB=((bm>>tv)&1)?-B[tv]:B[tv];
    }
    printf("NOFLIP minus=");
    for(int i=0;i<nb;++i)if((bm>>i)&1)printf("%d ",B[i]);
    printf("sigma_p=%+d TopPair=(%+d,%+d) parent=%s child=%s ratio=%s/%s\n",
      sp*p,((bm>>tu)&1)?-B[tu]:B[tu],((bm>>tv)&1)?-B[tv]:B[tv],
      pg.get_str().c_str(),cg.get_str().c_str(),
      cg.get_str().c_str(),pg.get_str().c_str());
  }
  printf("patterns=%lld noflip=%lld TopPair_fail=%lld max_ratio=%s/%s "
         "parent=%s child=%s TopPair=(%+d,%+d)\n",
      patterns,noflip,tpfail,maxnum.get_str().c_str(),maxden.get_str().c_str(),
      maxparent.get_str().c_str(),maxchild.get_str().c_str(),maxA,maxB);
}
'''

ENV=os.environ.copy()
ENV["TMPDIR"]="/dev/shm"
ENV["OMP_NUM_THREADS"]="4"

def compile_mem(name, source, omp=False):
    fd=os.memfd_create(name)
    os.set_inheritable(fd,True)
    cmd=["g++","-pipe","-std=c++17","-O3"]
    if omp: cmd.append("-fopenmp")
    cmd += ["-x","c++","-o",f"/proc/self/fd/{fd}","-","-lgmpxx","-lgmp"]
    q=subprocess.run(cmd,input=source,env=ENV,pass_fds=(fd,),
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if q.returncode: raise RuntimeError(q.stderr.decode())
    os.fchmod(fd,0o700)
    return fd

fast=compile_mem("fm174_fast",CPP.encode())

def scan_fast(profiles,name,progress_every):
    bands={b:{"profiles":0,"patterns":0,"noflip":0,"fail":0,
              "max":Fraction(-1),"record":None}
           for b in ("60-99","100-199","200-300")}
    recs=[]
    for i,(W,k,step,s,p,ns) in enumerate(profiles,1):
        if i==1 or i%progress_every==0 or i==len(profiles):
            print(time.strftime("%H:%M:%S"),name,"start",i,"/",len(profiles),
                  "W",W,"k",k,"p",p,flush=True)
        z=subprocess.run([f"/proc/self/fd/{fast}",str(p),
                          *[f"{n}:1" for n in ns]],
                         pass_fds=(fast,),env=ENV,text=True,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         timeout=90)
        if z.returncode: raise RuntimeError(z.stderr+z.stdout)
        band="60-99" if W<100 else ("100-199" if W<200 else "200-300")
        st=bands[band]; st["profiles"]+=1
        for line in z.stdout.splitlines():
            if line.startswith("NOFLIP minus="):
                m=re.search(r"minus=(.*?)sigma_p=([+-]?\d+) "
                            r"TopPair=\(([+-]?\d+),([+-]?\d+)\) "
                            r"parent=(\d+) child=(\d+)",line)
                if not m: raise ValueError(line)
                minus,sig,a,b,parent,child=m.groups()
                rec={"W":W,"k":k,"step":step,"s":s,"p":p,
                     "minus":minus.strip(),"sigma_p":int(sig),
                     "TopPair":(int(a),int(b)),"parent":int(parent),
                     "child":int(child),"ratio":Fraction(int(child),int(parent))}
                recs.append(rec);st["noflip"]+=1
                if rec["child"]>rec["parent"]:st["fail"]+=1
                if rec["ratio"]>st["max"]:
                    st["max"]=rec["ratio"];st["record"]=rec
                print("NOFLIP_RECORD",rec,flush=True)
            elif line.startswith("patterns="):
                m=re.search(r"patterns=(\d+) noflip=(\d+) TopPair_fail=(\d+)",line)
                st["patterns"]+=int(m.group(1));st["fail"]+=int(m.group(3))
        if i==len(profiles):
            print(time.strftime("%H:%M:%S"),name,"done",i,flush=True)
    print(name,"BANDS",bands)
    print(name,"TOTAL patterns",sum(v["patterns"] for v in bands.values()),
          "no-flip",len(recs),"TopPair failures",
          sum(v["fail"] for v in bands.values()))
    for v in bands.values():
        if v["record"] is not None:
            assert 16*v["max"].numerator < v["max"].denominator
    return bands,recs

canonical=[]
for k in range(8,17):
    for step in (1,2):
        for s in range(1,201):
            ns=[s+step*i for i in range(k)]
            W=sum(ns); mx=ns[-1]
            p=mx+2
            if (W-p)%2:p+=1
            delta=(W-p)//2
            if 60<=W<=300 and p>=6 and delta>=8 and mx<=delta \
                    and sum(x>=3 for x in ns)>=2:
                canonical.append((W,k,step,s,p,ns))
canonical.sort()
assert len(canonical)==313
canon_stats,canon_recs=scan_fast(canonical,"canonical-p",1)
assert [(canon_stats[b]["profiles"],canon_stats[b]["patterns"],
         canon_stats[b]["noflip"],canon_stats[b]["fail"])
        for b in ("60-99","100-199","200-300")] == [
            (29,24348,88,0),(127,308430,319,0),
            (157,536433,420,0)]

pmax_profiles=[]
for k in range(8,13):
    for step in (1,2):
        for s in range(1,201):
            ns=[s+step*i for i in range(k)]
            W=sum(ns);mx=ns[-1];need=max(8,mx)
            plo=max(6,mx)
            if (W-plo)%2:plo+=1
            p=W-2*need
            if (W-p)%2:p-=1
            if 60<=W<=300 and p>=plo and sum(x>=3 for x in ns)>=2:
                pmax_profiles.append((W,k,step,s,p,ns))
pmax_profiles.sort()
assert len(pmax_profiles)==230
pmax_stats,pmax_recs=scan_fast(pmax_profiles,"max-p",20)
assert [(pmax_stats[b]["profiles"],pmax_stats[b]["patterns"],
         pmax_stats[b]["noflip"],pmax_stats[b]["fail"])
        for b in ("60-99","100-199","200-300")] == [
            (28,20266,52,0),(98,92738,74,0),
            (104,102778,14,0)]

# Exact repository tools for the repeated-label profiles and independent checks.
root="ginibre_q3/character_ring_iter/fm39/"
def compile_repo(name,path,omp=False):
    return compile_mem(name,open(root+path,"rb").read(),omp)
many=compile_repo("fm174_many","sec166_class_pattern_flips.cpp",True)
flip=compile_repo("fm174_flip","sec166_flip_single.cpp",True)
gp=compile_repo("fm174_gp","sec162_gp_single.cpp",False)

seeds=[
 ("W48",{1:2,2:1,3:4,4:3,5:4},8),
 ("W50",{1:2,2:2,3:4,4:3,5:4},6),
 ("W56",{1:2,2:1,3:4,4:3,5:4,8:1},8)]

def pgrid(d):
    W=sum(n*m for n,m in d.items());mx=max(d);need=max(8,mx)
    lo=max(6,mx)
    if (W-lo)%2:lo+=1
    hi=W-2*need
    if (W-hi)%2:hi-=1
    if hi<lo:return []
    mid=lo+2*((hi-lo)//4)
    return sorted(set((lo,mid,hi)))

profiles={}
def add_profile(tag,d,p):
    d={int(n):int(m) for n,m in d.items() if m>0}
    W=sum(n*m for n,m in d.items());mx=max(d)
    delta=(W-p)//2
    if 48<=W<=300 and p>=max(6,mx) and (W-p)%2==0 \
       and delta>=8 and mx<=delta \
       and sum(m for n,m in d.items() if n>=3)>=2:
        items=tuple(sorted(d.items()))
        profiles[(items,p)]=tag

for name,base,special in seeds:
    for t in range(1,6):
        d={n:m*t for n,m in base.items()}
        ps=pgrid(d) if t<=3 else [max(pgrid(d))]
        for p in ps:add_profile(f"{name}*{t}",d,p)
        if t==1:add_profile(f"{name}-known",d,special)

base48=seeds[0][1]
for t,extras in [(1,(6,)),(1,(7,)),(1,(8,)),(1,(6,7,8)),
                 (2,(6,)),(2,(7,)),(2,(8,))]:
    d={n:m*t for n,m in base48.items()}
    for n in extras:d[n]=d.get(n,0)+1
    ps=pgrid(d)
    if ps:add_profile("W48*%d+%s"%(t,"-".join(map(str,extras))),
                      d,ps[len(ps)//2])

ordered=sorted((sum(n*m for n,m in items),p,items,tag)
               for (items,p),tag in profiles.items())
assert len(ordered)==40

def candidate_pattern_count(items,p):
    d=dict(items);labels=list(d);mults=list(d.values());out=0
    for mask in range(1<<len(labels)):
        neg=sum(mults[i] for i in range(len(labels)) if mask>>i&1)
        if not 2<=neg<=6:continue
        sp=-1 if neg%2 else 1
        if any(labels[i]==p and
               (-1 if mask>>i&1 else 1)!=sp
               for i in range(len(labels))):continue
        out+=1
    return out

shape_bands={b:{"runs":0,"patterns":0,"noflip":0,"fail":0,
                "max":Fraction(-1),"record":None}
             for b in ("48-99","100-199","200-300")}
shape_records=[]
for i,(W,p,items,tag) in enumerate(ordered,1):
    d=dict(items);band="48-99" if W<100 else ("100-199" if W<200 else "200-300")
    st=shape_bands[band];st["runs"]+=1
    st["patterns"]+=candidate_pattern_count(items,p)
    print(time.strftime("%H:%M:%S"),"shape start",i,"/",len(ordered),
          tag,"W",W,"p",p,flush=True)
    z=subprocess.run([f"/proc/self/fd/{many}",str(p),
                      *[f"{n}:{m}" for n,m in items]],
                     pass_fds=(many,),env=ENV,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
    if z.returncode:raise RuntimeError(z.stderr+z.stdout)
    for line in z.stdout.splitlines():
        if not line.startswith("NOFLIP pattern"):continue
        m=re.search(r"minus-classes=(.*?) sigma p=([+-]?\d+) "
                    r"TopPair\(([-+]?\d+),([-+]?\d+)\) "
                    r"(ok|FAIL) g=(\d+) child=(\d+)",line)
        if not m:raise ValueError(line)
        minus,sig,a,b,status,parent,child=m.groups()
        neg=sum(d[int(x)] for x in minus.split()) if minus.strip() else 0
        if not 2<=neg<=6:continue
        rec={"tag":tag,"W":W,"p":p,"delta":(W-p)//2,
             "items":items,"minus":minus.strip(),"sigma_p":int(sig),
             "TopPair":(int(a),int(b)),"parent":int(parent),
             "child":int(child),"ratio":Fraction(int(child),int(parent))}
        shape_records.append(rec);st["noflip"]+=1
        if status=="FAIL":st["fail"]+=1
        if rec["ratio"]>st["max"]:st["max"]=rec["ratio"];st["record"]=rec
        print("SHAPE_NOFLIP_RECORD",rec,flush=True)
    print(time.strftime("%H:%M:%S"),"shape done",i,flush=True)

print("SHAPE_BANDS",shape_bands)
print("SHAPE candidate patterns",
      sum(v["patterns"] for v in shape_bands.values()),
      "no-flip",len(shape_records),
      "TopPair failures",sum(v["fail"] for v in shape_bands.values()))
assert sum(v["patterns"] for v in shape_bands.values())==339
assert len(shape_records)==3
assert all(r["child"]<=r["parent"] for r in shape_records)
assert max(r["ratio"] for r in shape_records)==Fraction(84983,2125083)
assert 16*84983<2125083

def check_repo(B,p,rem,parent,child):
    sigma=-1 if sum(x<0 for x in B)%2 else 1
    z=subprocess.run([f"/proc/self/fd/{flip}",str(sigma*p),*map(str,B)],
                     pass_fds=(flip,),env=ENV,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
    assert z.returncode==0 and "NO flip descent among" in z.stdout
    q=subprocess.run([f"/proc/self/fd/{gp}",str(p),*map(str,B),"--",*map(str,rem)],
                     pass_fds=(gp,),env=ENV,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
    a=re.search(r"g_p\(B\) = (\d+)",q.stdout)
    b=re.search(r"g_p\(B-R\) = (\d+)",q.stdout)
    assert q.returncode==0 and (int(a.group(1)),int(b.group(1)))==(parent,child)
    print("REPO_CROSSCHECK",p,B,rem,parent,child)

checks=[
 (13,[3,-4,5,-6,7,-8,9,-10,11],[9,11],328561,4444),
 (16,[2,-3,4,5,6,-7,8,9,10,11,12,13,14],[12,14],
  1305099613,10392078),
 (22,[-5,6,-7,8,9,10,11,12,13,14,15,16,17,18,19,20],[18,20],
  318194106750922,1131644431324),
 (38,[4,-5,-6,7,-8,9,-10,11],[9,11],10594,5),
 (64,[6,8,10,-12,-14,-16,18,20],[18,20],183406,5),
 (166,[14,16,-18,20,-22,-24,-26,-28,-30,32],[-30,32],75019832,7),
 (8,[1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5],[5,5],4075371,160741),
 (6,[1,1,-2,-2,3,3,3,3,-4,-4,-4,5,5,5,5],[5,5],10625415,424915),
 (8,[1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8],[-4,8],27576641,961895)]
for case in checks:check_repo(*case)
print("PASS")
