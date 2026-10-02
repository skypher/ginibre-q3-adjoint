import ctypes, os, pathlib, random, re, subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from math import comb

def sparse(B, tilted):
    P={(0,0,0):1}
    for z in B:
        n=abs(z); e=1 if z>0 else -1; T=defaultdict(int)
        for j in range(n+1): T[j,j,0]+=1
        if tilted:
            for r in range(n//2+1):
                h=n-2*r
                for j in range(h+1):
                    T[j+r,h-j+r,h]+=e*(-1)**r*comb(n-r,r)*comb(h,j)
        else:
            for j in range(n+1): T[j,n-j,0]+=e
        N=defaultdict(int)
        for (i,j,k),v in P.items():
            for (a,b,c),w in T.items(): N[i+a,j+b,k+c]+=v*w
        P={x:y for x,y in N.items() if y}
    return P

def hval(M,i,j):
    if i<0 or j<0: return 0
    if i<j: i,j=j,i
    return M.get((i,j,0),0)-M.get((i+1,j-1,0),0)

def jrec(M,k,i,j):
    if i<0 or j<0: return Fraction(0)
    if i<j: i,j=j,i
    alpha=Fraction(k-i+j,2); beta=Fraction(k+i-j,2)
    c=1/(beta+1); out=Fraction(0)
    for t in range(j+1):
        out+=c*hval(M,i+t,j-t)
        c*=(alpha-t)/(beta+t+2)
    return out

moment_checks=0
for B in ((1,-2,3),(-1,2,-3,4),(1,1,-2,3,-4),(-2,-4,5),(-1,1,2,-2)):
    M=sparse(B,False); V=sparse(B,True); W=sum(map(abs,B))
    for k in range(8):
        for i in range(W+1):
            for j in range(i+1):
                direct=sum((Fraction(2*v,k+r+2)
                            for (x,y,r),v in V.items() if (x,y)==(i,j)),
                           Fraction(0))
                assert jrec(M,k,i,j)==direct
                moment_checks+=1
assert moment_checks==2224
print("direct moment recurrence checks:",moment_checks,"PASS")

base=pathlib.Path(
    "/tmp/claude-1006/-home-yang-q3adjoint/"
    "2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/st52")
pat=re.compile(
    r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+) D=(\w+)"
    r" TopPair\((-?\d+),(-?\d+)\)=(\w+)")
rows=[]; byW=Counter()
for done in sorted(base.glob("chunk_*.done"),
                   key=lambda f:int(f.stem.split("_")[1])):
    k=int(done.stem.split("_")[1]); src=base/f"chunk_{k}.noflip"; part=[]
    if src.exists():
        for line in src.read_text().splitlines():
            m=pat.search(line)
            if not m: raise ValueError(("line parse failed",k,line))
            w,ps,bs,phi,dstat,ta,tb,tp=m.groups()
            B=tuple(map(int,bs.split()))
            if dstat!="ok" or sum(map(abs,B))!=int(w):
                raise ValueError(("bad row",line))
            part.append((abs(int(ps)),int(phi),int(ta),int(tb),
                         len(B),int(ps),B))
            byW[int(w)]+=1
    rows.extend(part)
    print(f"[{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S UTC}] "
          f"chunk_{k}: {len(part)} rows; cumulative={len(rows)}",flush=True)
temps=sorted(x.name for x in base.glob("chunk_*.tmp"))
print("temporary chunks:",temps)
assert rows

rng=random.Random(104)
picks=sorted(rng.sample(range(len(rows)),min(512,len(rows))))
items=[(0,*r) for r in rows]+[(1,*rows[i]) for i in picks]
data=("\n".join([str(len(items))]+[
    f"{tag} {p} {phi} {ta} {tb} {n} {ps} "+" ".join(map(str,B))
    for tag,p,phi,ta,tb,n,ps,B in items])+"\n").encode()

CPP=r'''
#include <bits/stdc++.h>
#include <gmpxx.h>
#include <omp.h>
using namespace std; using I=mpz_class; using Q=mpq_class;
struct Row{int tag,p,phi,ta,tb,n,ps;vector<int>B;};
struct Tab{int D;vector<I>v;I at(int i,int j)const{
 return i<0||j<0||i>=D||j>=D?I(0):v[(size_t)i*D+j];}};
int weight(const vector<int>&B){int s=0;for(int z:B)s+=abs(z);return s;}

Tab fuse(const vector<int>&B){
 int D=weight(B)+1;size_t z=(size_t)D*D;
 vector<I>f(z),g(z);vector<int>active{0};f[0]=1;
 for(int x:B){
  int n=abs(x),e=x>0?1:-1;vector<int>touch;vector<unsigned char>seen(z);
  for(int id:active){
   int i=id/D,j=id%D;const I&v=f[id];if(v==0)continue;
   for(int a=abs(i-n);a<=i+n;a+=2){
    int q=a*D+j;if(!seen[q])seen[q]=1,touch.push_back(q);g[q]+=v;
   }
   for(int b=abs(j-n);b<=j+n;b+=2){
    int q=i*D+b;if(!seen[q])seen[q]=1,touch.push_back(q);g[q]+=e*v;
   }
  }
  vector<int>next;for(int q:touch)if(g[q]!=0)next.push_back(q);
  for(int q:active)f[q]=0;
  f.swap(g);active.swap(next);
 }
 return {D,std::move(f)};
}

Tab ordinary(const vector<int>&B){
 int D=weight(B)+1;size_t z=(size_t)D*D;
 vector<I>f(z),g(z);vector<int>active{0};f[0]=1;
 for(int x:B){
  int n=abs(x),e=x>0?1:-1;vector<int>touch;vector<unsigned char>seen(z);
  for(int id:active){
   int i=id/D,j=id%D;const I&v=f[id];if(v==0)continue;
   for(int t=0;t<=n;t++){
    int q=(i+t)*D+j+t;
    if(!seen[q])seen[q]=1,touch.push_back(q);g[q]+=v;
    q=(i+t)*D+j+n-t;
    if(!seen[q])seen[q]=1,touch.push_back(q);g[q]+=e*v;
   }
  }
  vector<int>next;for(int q:touch)if(g[q]!=0)next.push_back(q);
  for(int q:active)f[q]=0;
  f.swap(g);active.swap(next);
 }
 return {D,std::move(f)};
}

I ceilroot(const Q&x){
 if(x<0)throw runtime_error("negative E");
 I n=x.get_num(),d=x.get_den(),u=n/d,r;
 mpz_sqrt(r.get_mpz_t(),u.get_mpz_t());
 if(r*r*d<n)++r;
 return r;
}
struct Cert{bool pass;I parent,child,margin;Q q,E,ratio;int k;};

Cert certificate(const vector<int>&B,int p,const Tab&FB,int ia,int ib){
 int delta=(weight(B)-p)/2,za=B[ia],zb=B[ib];
 if(abs(za)<abs(zb))swap(za,zb);
 int a=abs(za),b=abs(zb),ea=za>0?1:-1,eb=zb>0?1:-1;
 vector<int>C;
 for(int i=0;i<(int)B.size();i++)if(i!=ia&&i!=ib)C.push_back(B[i]);
 Tab F=fuse(C),M=ordinary(C);
 auto H=[&](int i,int j)->I{
  if(i<0||j<0)return 0;
  if(i<j)swap(i,j);
  return M.at(i,j)-M.at(i+1,j-1);
 };
 map<tuple<int,int,int>,Q>cache;
 function<Q(int,int,int)>J=[&](int k,int i,int j)->Q{
  if(i<0||j<0)return 0;
  if(i<j)swap(i,j);
  auto key=make_tuple(k,i,j);auto at=cache.find(key);
  if(at!=cache.end())return at->second;
  Q alpha(k-i+j,2),beta(k+i-j,2);
  alpha.canonicalize();beta.canonicalize();
  Q c=Q(1)/(beta+1),v=0;
  for(int t=0;t<=j;t++){
   v+=c*H(i+t,j-t);
   c*=(alpha-t)/(beta+t+2);
  }
  cache.emplace(key,v);return v;
 };
 auto g=[&](int i,int j)->I{return F.at(i,j);};
 I S=0,Z=0;
 for(int c=a-b;c<=a+b;c+=2){
  for(int u=abs(p-c);u<=p+c;u+=2)S+=g(u,0);
  Z+=ea*eb*g(p,c);
 }
 int r=delta-a,s=delta-b,ell=delta-a-b-1;
 I A=eb*H(delta,s)+ea*H(delta,r);
 I wrap=eb*H(r-1,ell)+ea*H(s-1,ell);
 I q=S+Z-g(p,0)-wrap,parent=FB.at(p,0),child=g(p,0);
 if(parent-child!=q+A)throw runtime_error("parent-child identity");
 Q best;bool first=true;int bestk=-1;
 for(int k=0;k<=2*b;k++){
  Q V=Q((b+1)*(b+1))*J(2*b-k,s,s)
     +Q((a+1)*(a+1))*J(2*a-k,r,r)
     +Q(2*ea*eb*(a+1)*(b+1))*J(a+b-k,r,s);
  Q jd=J(k,delta,delta),E=jd*V;
  if(jd<0||V<0||E<0||Q(A*A)>E)
   throw runtime_error("Gram/Cauchy check");
  if(first||E<best){best=E;bestk=k;first=false;}
 }
 bool pass=q>=0&&Q(q*q)>=best;
 I margin=q-ceilroot(best);
 Q ratio=parent>0?Q(child,parent):Q(0);
 return {pass,parent,child,margin,Q(q),best,ratio,bestk};
}

string qs(const Q&x){return x.get_str();}
string is(const I&x){return x.get_str();}
string vecstr(const vector<int>&v){
 ostringstream o;o<<"(";
 for(size_t i=0;i<v.size();i++){if(i)o<<",";o<<v[i];}
 o<<")";return o.str();
}

extern "C" int audit(const char* bytes){
 istringstream in(bytes);int n;
 if(!(in>>n)||n<0||n>200000)return 2;
 vector<Row>R(n);
 for(auto&r:R){
  in>>r.tag>>r.p>>r.phi>>r.ta>>r.tb>>r.n>>r.ps;
  if(!in||r.n<2||r.n>100)return 3;
  r.B.resize(r.n);
  for(int&x:r.B){in>>x;if(!in||x==0)return 4;}
 }
 vector<Cert>ans(n);vector<string>fail(n);
 atomic<int>errors{0},done{0},samples{0};
 #pragma omp parallel for schedule(dynamic,4)
 for(int ix=0;ix<n;ix++){
  Row&r=R[ix];
  try{
   if(r.tag==1){
    int neg=0;for(int x:r.B)neg+=x<0;
    int sig=neg%2?-1:1;
    if(r.ps!=sig*r.p)throw runtime_error("suffix parity");
    vector<int>L=r.B;L.push_back(sig*r.p);bool descent=false;
    for(int i=0;i<(int)L.size();i++)for(int j=i+1;j<(int)L.size();j++){
     vector<int>C;
     for(int t=0;t<(int)L.size();t++)if(t!=i&&t!=j)C.push_back(L[t]);
     Tab T=fuse(C);
     I x=(L[j]>0?1:-1)*T.at(abs(L[i]),abs(L[j]));
     I y=(L[i]>0?1:-1)*T.at(abs(L[j]),abs(L[i]));
     if(x!=y)throw runtime_error("flip orientation mismatch");
     if(x>=0)descent=true;
    }
    if(descent)throw runtime_error("sample not no-flip");
    ++samples;continue;
   }
   Tab FB=fuse(r.B);
   if(2*FB.at(r.p,0)!=r.phi)throw runtime_error("phi mismatch");
   int ia=-1,ib=-1,bw=-1,bm=-1;
   for(int i=0;i<(int)r.B.size();i++)for(int j=i+1;j<(int)r.B.size();j++)
    if((abs(r.B[i])+abs(r.B[j]))%2==0){
     int w=abs(r.B[i])+abs(r.B[j]),m=max(abs(r.B[i]),abs(r.B[j]));
     if(w>bw||(w==bw&&m>bm)){ia=i;ib=j;bw=w;bm=m;}
    }
   if(ia<0||!((r.B[ia]==r.ta&&r.B[ib]==r.tb)||
              (r.B[ia]==r.tb&&r.B[ib]==r.ta))
      )throw runtime_error("TopPair mismatch");
   Cert top=certificate(r.B,r.p,FB,ia,ib);ans[ix]=top;
   if(!top.pass){
    ostringstream o;
    o<<"W="<<weight(r.B)<<" p="<<r.ps<<" B="<<vecstr(r.B)
     <<" TP=("<<r.B[ia]<<","<<r.B[ib]<<") g="<<is(top.parent)<<"/"
     <<is(top.child)<<" ratio="<<qs(top.ratio)<<" Q="<<qs(top.q)
     <<" E_min="<<qs(top.E)<<" k="<<top.k
     <<" margin="<<is(top.margin)<<"\n";
    bool any=false;int nr=0;
    for(int i=0;i<(int)r.B.size();i++)for(int j=i+1;j<(int)r.B.size();j++)
     if((abs(r.B[i])+abs(r.B[j]))%2==0){
      ++nr;Cert c=certificate(r.B,r.p,FB,i,j);if(c.pass)any=true;
      o<<" remove=("<<r.B[i]<<","<<r.B[j]<<") Q="<<qs(c.q)
       <<" E_min="<<qs(c.E)<<" k="<<c.k<<" margin="<<is(c.margin)
       <<" pass="<<c.pass<<"\n";
     }
    o<<"same_parity_removals="<<nr<<" any_pass="<<any;
    if(!any)o<<" ALL_REMOVALS_FAIL";
    fail[ix]=o.str();
   }
   int d=++done;
   if(d%10000==0){
    #pragma omp critical
    {
     time_t z=time(nullptr);tm*t=gmtime(&z);char b[40];
     strftime(b,sizeof(b),"%Y-%m-%d %H:%M:%S UTC",t);
     cerr<<"["<<b<<"] checked "<<d<<" rows\n";
    }
   }
  }catch(const exception&e){
   ++errors;fail[ix]=string("ERROR row ")+to_string(ix)+": "+e.what();
  }
 }
 int total=0,passes=0,mono=0,nonmono=0,allfail=0,rescued=0,wi=-1,ri=-1;
 I worstmargin;Q maxratio;bool have=false;
 for(int i=0;i<n;i++)if(R[i].tag==0){
  ++total;Cert&c=ans[i];passes+=c.pass;
  if(c.parent>=c.child)++mono;else ++nonmono;
  if(!have||c.ratio>maxratio){maxratio=c.ratio;ri=i;have=true;}
  if(wi<0||c.margin<worstmargin){wi=i;worstmargin=c.margin;}
  if(!fail[i].empty()&&fail[i].find("ERROR ")!=0){
   if(fail[i].find("any_pass=0")!=string::npos)++allfail;else ++rescued;
  }
 }
 auto topidx=[&](const vector<int>&B){
  int x=-1,y=-1,bw=-1,bm=-1;
  for(int i=0;i<(int)B.size();i++)for(int j=i+1;j<(int)B.size();j++)
   if((abs(B[i])+abs(B[j]))%2==0){
    int w=abs(B[i])+abs(B[j]),m=max(abs(B[i]),abs(B[j]));
    if(w>bw||(w==bw&&m>bm)){x=i;y=j;bw=w;bm=m;}
   }
  return pair<int,int>(x,y);
 };
 time_t z=time(nullptr);tm*t=gmtime(&z);char stamp[40];
 strftime(stamp,sizeof(stamp),"%Y-%m-%d %H:%M:%S UTC",t);
 cout<<"["<<stamp<<"] records="<<total
     <<" TopPair_CB_pass="<<passes<<" fail="<<total-passes
     <<" alt_rescued="<<rescued<<" fail_all_removals="<<allfail<<"\n";
 cout<<"TopPair_monotone="<<mono<<" nonmonotone="<<nonmono
     <<" worst_ratio="<<qs(maxratio)<<"\n";
 for(int ix:{ri,wi})if(ix>=0){
  auto [a,b]=topidx(R[ix].B);Cert&c=ans[ix];
  cout<<(ix==ri?"max_ratio_witness ":"min_margin_witness ")
      <<"W="<<weight(R[ix].B)<<" p="<<R[ix].ps
      <<" B="<<vecstr(R[ix].B)<<" TP=("<<R[ix].B[a]<<","<<R[ix].B[b]<<")"
      <<" g="<<is(c.parent)<<"/"<<is(c.child)<<" ratio="<<qs(c.ratio)
      <<" Q="<<qs(c.q)<<" E_min="<<qs(c.E)<<" k="<<c.k
      <<" margin="<<is(c.margin)<<"\n";
 }
 cout<<"random_no_flip="<<samples<<"/512 seed=104\n";
 for(auto&s:fail)if(!s.empty())cout<<s<<"\n";
 cout<<"errors="<<errors<<"\n";
 return errors?1:0;
}
'''
obj=os.memfd_create("fm104-verified-object")
lib=os.memfd_create("fm104-verified-library")
subprocess.run(["g++","-std=c++17","-O2","-fPIC","-fopenmp","-pipe","-x","c++","-c","-",
                "-o",f"/proc/self/fd/{obj}"],
               input=CPP,text=True,pass_fds=(obj,),check=True)
subprocess.run(["ld","-shared",f"/proc/self/fd/{obj}",
                "/lib/x86_64-linux-gnu/libgmpxx.so.4",
                "/lib/x86_64-linux-gnu/libgmp.so.10",
                "/lib/x86_64-linux-gnu/libstdc++.so.6",
                "/lib/x86_64-linux-gnu/libgomp.so.1","-lc",
                "-o",f"/proc/self/fd/{lib}"],pass_fds=(obj,lib),check=True)
dll=ctypes.CDLL(f"/proc/self/fd/{lib}")
dll.audit.argtypes=[ctypes.c_char_p];dll.audit.restype=ctypes.c_int
assert dll.audit(data)==0
os.close(obj);os.close(lib)
print("rows by weight bands:",
      {f"{x}-{x+9}":sum(v for w,v in byW.items() if x<=w<=x+9)
       for x in (23,33,43)})
