import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import random

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

def entry(B,p,q=0):
    if min(p,q)<0 or p+q>sum(map(abs,B)): return 0
    states={(0,0):1}
    rem=sum(map(abs,B))
    for z in sorted(B,key=abs,reverse=True):
        n=abs(z); eps=1 if z>0 else -1; rem-=n
        out=defaultdict(int)
        for (s,t),v in states.items():
            r=rem-abs(t-q)
            if r>=0:
                lo=max(abs(s-n),p-r)
                lo+=(s+n-lo)%2
                for h in range(lo,min(s+n,p+r)+1,2):
                    out[h,t]+=v
            r=rem-abs(s-p)
            if r>=0:
                lo=max(abs(t-n),q-r)
                lo+=(t+n-lo)%2
                for h in range(lo,min(t+n,q+r)+1,2):
                    out[s,h]+=eps*v
        states={st:v for st,v in out.items() if v}
    return states.get((p,q),0)

def table(B):
    states={(0,0):1}
    for z in B:
        out=defaultdict(int);n=abs(z);eps=1 if z>0 else -1
        for (s,t),v in states.items():
            for h in cg(s,n): out[h,t]+=v
            for h in cg(t,n): out[s,h]+=eps*v
        states={st:v for st,v in out.items() if v}
    return states

def top(B):
    i,j=max(((i,j) for i in range(len(B)) for j in range(i+1,len(B))
             if (B[i]-B[j])%2==0),
            key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                           max(abs(B[ij[0]]),abs(B[ij[1]]))))
    return (B[i],B[j]),tuple(z for k,z in enumerate(B) if k not in (i,j))

def budget(B,p,ref=None):
    counts=Counter(map(abs,B));refcounts=Counter(map(abs,B if ref is None else ref))
    ns=tuple(sorted(refcounts));ms=tuple(counts[n] for n in ns)
    @lru_cache(None)
    def mult(ct):
        if not any(ct): return (1,)
        i=next(i for i,m in enumerate(ct) if m)
        prev=list(ct);prev[i]-=1;ns_i=ns[i]
        v=mult(tuple(prev));out=[0]*(len(v)+ns_i)
        for s,x in enumerate(v):
            if x:
                for h in cg(s,ns_i): out[h]+=x
        return tuple(out)
    total=forced=0
    parity=sum((refcounts[n]%2)<<i for i,n in enumerate(ns))
    constraint=(parity^(1<<ns.index(p))) if p in ns else None
    H=[0]*(len(B)+1)
    for ct in product(*(range(m+1) for m in ms)):
        co=tuple(m-c for m,c in zip(ms,ct));v=mult(co)
        if p>=len(v): continue
        choose=1;mask=0
        for i,(m,c) in enumerate(zip(ms,ct)):
            choose*=comb(m,c);mask|=(c%2)<<i
        w=choose*mult(ct)[0]*v[p]
        total+=w;H[sum(ct)]+=w
        if mask==0 or mask==constraint: forced+=w
    return total,forced,H

def pieces(B,p):
    R,C=top(B);a,b=map(abs,R);ea=1 if R[0]>0 else -1;eb=1 if R[1]>0 else -1
    T=table(C);get=lambda s,t:T.get((s,t),0)
    same=sum(get(s,0) for c in cg(a,b) for s in cg(p,c))
    x=eb*sum(get(s,b) for s in cg(p,a))
    y=ea*sum(get(s,a) for s in cg(p,b))
    z=ea*eb*sum(get(p,c) for c in cg(a,b))
    return same,x,y,z,get(p,0)

def edge(C,d,eps):
    p=sum(map(abs,C));assert d<=p
    T=table(C)
    P=[sum(v for (s,t),v in T.items()
           if t==0 and abs(p-2*j)<=s and (s-p)%2==0) for j in range(d+1)]
    E=[1]+[0]*d
    for z in C:
        n=abs(z);sg=1 if z>0 else -1
        for j in range(d,n-1,-1): E[j]+=sg*E[j-n]
    e=abs(E[d]);L=d+1
    assert all(x>=1 for x in P)
    assert L*P[d]>e*e
    g=entry(tuple(C)+(eps*d,eps*d),p)
    assert g==1+sum(P)+2*eps*E[d]
    assert g-1 >= (e-L)**2//L
    return g,E[d],P

def full_pairfree(B,p):
    sig=(-1)**sum(z<0 for z in B)
    L=tuple(B)+(sig*p,)
    return all(-z not in L for z in L)

def noflip(B,p):
    sig=(-1)**sum(z<0 for z in B);L=tuple(B)+(sig*p,)
    ans=[]
    for i in range(len(L)):
        for j in range(i+1,len(L)):
            C=L[:i]+L[i+1:j]+L[j+1:]
            v=(1 if L[j]>0 else -1)*entry(C,abs(L[i]),abs(L[j]))
            ans.append(v)
    return max(ans),ans

CPP=r'''
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using I=__int128_t;
I pw(I a,int n){I z=1;while(n-->0)z*=a;return z;}
I choose(int n,int k){I z=1;for(int j=1;j<=k;j++)z=z*(n-j+1)/j;return z;}
vector<I> table(const vector<int>&B){
 int W=0;for(int z:B)W+=abs(z);int D=W+1,s=0;
 vector<I>A(D*D),Q(D*D);A[0]=1;
 for(int z:B){int n=abs(z),e=z>0?1:-1;fill(Q.begin(),Q.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){I v=A[i*D+j];if(!v)continue;
   for(int h=abs(i-n);h<=i+n;h+=2)Q[h*D+j]+=v;
   for(int h=abs(j-n);h<=j+n;h+=2)Q[i*D+h]+=e*v;}
  A.swap(Q);s+=n;}return A;
}
struct Masks{vector<I>B,C,H;int forced=-1;};
Masks masks(const vector<int>&B,const vector<int>&C,int p){
 map<int,int>bm,cm;for(int z:B)bm[abs(z)]++;for(int z:C)cm[abs(z)]++;
 vector<int>ns,ms,cs,st;int L=1,W=0,ccode=0,vm=0,pk=-1;
 for(auto [n,m]:bm){int i=ns.size();ns.push_back(n);ms.push_back(m);
  cs.push_back(cm[n]);st.push_back(L);ccode+=cm[n]*L;L*=m+1;W+=n*m;
  vm^=(m%2)<<i;if(n==p)pk=i;}
 int K=ns.size(),D=W+1;
 vector<I>U((long long)L*D);U[0]=1;vector<int>wt(L);
 for(int code=1;code<L;code++){int i=0;while(code/st[i]%(ms[i]+1)==0)++i;
  int prev=code-st[i],n=ns[i];wt[code]=wt[prev]+n;
  for(int j=0;j<=wt[prev];j++){I v=U[(long long)prev*D+j];if(!v)continue;
   for(int h=abs(j-n);h<=j+n;h+=2)U[(long long)code*D+h]+=v;}}
 Masks out;out.B.resize(1<<K);out.C.resize(1<<K);out.H.resize(B.size()+1);
 if(pk>=0)out.forced=vm^(1<<pk);
 for(int code=0;code<L;code++){I cb=1,cc=1;int mask=0,size=0;
  for(int i=0;i<K;i++){int a=code/st[i]%(ms[i]+1);size+=a;mask|=(a%2)<<i;
   cb*=choose(ms[i],a);if(a>cs[i])cc=0;else cc*=choose(cs[i],a);}
  I b=p<=W?cb*U[(long long)code*D]*U[(long long)(L-1-code)*D+p]:0;
  out.B[mask]+=b;out.H[size]+=b;
  if(cc&&p<=W)out.C[mask]+=cc*U[(long long)code*D]*U[(long long)(ccode-code)*D+p];
 }return out;
}
extern "C" int audit(const char*data){
 istringstream in(data);int K;in>>K;vector<vector<int>>bs(K);vector<int>ps(K),tag(K);
 for(int i=0;i<K;i++){int n;in>>tag[i]>>ps[i]>>n;bs[i].resize(n);for(int&z:bs[i])in>>z;}
 array<atomic<long long>,2>cases{},old{},parity{},edge{},tail{},covered{},absolute{};
 atomic<long long>bad{0};
 #pragma omp parallel for schedule(dynamic,8)
 for(int q=0;q<K;q++){
  auto B=bs[q];int p=ps[q],N=B.size(),W=0,M=0,m=1000;
  for(int z:B){W+=abs(z);M=max(M,abs(z));m=min(m,abs(z));}
  if(W>40||N<5){++bad;continue;}
  int i0=-1,j0=-1,w=-1,ma=-1;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)if((B[i]-B[j])%2==0){
   int v=abs(B[i])+abs(B[j]),a=max(abs(B[i]),abs(B[j]));
   if(make_pair(v,a)>make_pair(w,ma))w=v,ma=a,i0=i,j0=j;}
  auto C=B;C.erase(C.begin()+j0);C.erase(C.begin()+i0);
  int a=abs(B[i0]),b=abs(B[j0]),ea=B[i0]>0?1:-1,eb=B[j0]>0?1:-1;
  int T=W-w,D=T+1,d=(W-p)/2;auto G=table(C),GP=table(B);
  auto at=[&](int s,int t)->I{return s+t<=T?G[s*D+t]:0;};
  I child=at(p,0),parent=GP[p*(W+1)],S=0,X=0,Y=0,Z=0;
  for(int c=abs(a-b);c<=a+b;c+=2){
   for(int s=abs(p-c);s<=p+c;s+=2)S+=at(s,0);
   Z+=ea*eb*at(p,c);}
  for(int s=abs(p-a);s<=p+a;s+=2)X+=eb*at(s,b);
  for(int s=abs(p-b);s<=p+b;s+=2)Y+=ea*at(s,a);
  if(parent!=S+X+Y+Z||S<(min(a,b)+1)*child)++bad;
  for(int s=0;s<=T;s++)if(at(s,0)<0)++bad;
  auto V=masks(B,C,p);I AB=0,AC=0,FB=0,FC=0;
  for(int mask=0;mask<(int)V.B.size();mask++){
   if(V.B[mask]<(min(a,b)+1)*V.C[mask])++bad;
   AB+=V.B[mask];AC+=V.C[mask];
   if(mask==0||mask==V.forced)FB+=V.B[mask],FC+=V.C[mask];}
  I lower=2*(FB-FC)-(AB-AC);
  if(parent<2*FB-AB||child>AC||parent-child<lower)++bad;
  bool pp=lower>=0,oo=2*FB>=AB+AC,ee=w==2*d;
  int h=min(m/2,d/(N-1)),dim=M+1;
  I badbound=0,badsum=0;
  if(V.H[0]<pw(h+1,N-2))++bad;
  for(int j=3;j<=N-2;j++){
   I z=choose(N,j)*pw(dim,N-5);badbound+=z;badsum+=V.H[j];
   if(V.H[j]>z)++bad;}
  if(parent<V.H[0]-badsum||AC>pw(2,N-2)*pw(dim,N-4))++bad;
  bool tt=pw(h+1,N-2)>=badbound+pw(2,N-2)*pw(dim,N-4);
  if(ee){
   vector<I>P(d+1),E(d+1);E[0]=1;
   for(int j=0;j<=d;j++)for(int s=abs(p-2*j);s<=p;s+=2)P[j]+=at(s,0);
   for(int z:C)for(int j=d;j>=abs(z);j--)E[j]+=(z>0?1:-1)*E[j-abs(z)];
   I es=E[d],e=es<0?-es:es,sp=0,L=d+1;
   for(I v:P){sp+=v;if(v<1)++bad;}
   if(L*P[d]<=e*e||parent!=1+sp+2*ea*es||parent-1<(e-L)*(e-L)/L)++bad;
  }
  auto ab=[](I x){return x<0?-x:x;};
  int t=tag[q];++cases[t];if(oo)++old[t];if(pp)++parity[t];if(ee)++edge[t];
  if(tt)++tail[t];if(pp||ee||tt)++covered[t];if(S-child>=ab(X)+ab(Y)+ab(Z))++absolute[t];
  if((pp||ee||tt)&&parent<child)++bad;
 }
 for(int t=0;t<2;t++)cout<<"tag="<<t<<" cases="<<cases[t]<<" old="<<old[t]
  <<" parity="<<parity[t]<<" boundary="<<edge[t]<<" tail="<<tail[t]
  <<" union="<<covered[t]<<" abs_blocks="<<absolute[t]<<"\n";
 cout<<"audit errors="<<bad<<"\n"<<flush;return bad?1:0;
}
'''

def check_bounds(B,p):
    R,C=top(B);A,F,H=budget(B,p);AC,FC,HC=budget(C,p,B)
    g=entry(B,p);gc=entry(C,p)
    S,X,Y,Z,child=pieces(B,p)
    assert S+X+Y+Z==g and child==gc
    assert S>=(min(map(abs,R))+1)*gc
    assert g>=2*F-A and gc<=AC
    assert g-gc>=2*(F-FC)-(A-AC)
    N=len(B);M=max(map(abs,B));d=(sum(map(abs,B))-p)//2
    h=min(min(map(abs,B))//2,d//(N-1))
    assert H[0]>=(h+1)**(N-2)
    assert all(H[j]<=comb(N,j)*(M+1)**(N-5) for j in range(3,N-1))
    assert g>=H[0]-sum(H[3:N-1])
    assert AC<=2**(N-2)*(M+1)**(N-4)
    return g,gc,A,F,AC,FC

def census(path):
    import ctypes,os,subprocess
    rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('NOFLIP','TPFAIL')): continue
        p=abs(int(line.split(' p=')[1].split()[0]))
        B=list(map(int,line.split(' B=')[1].split(' phi=')[0].split(' flip=')[0].split()))
        rows.append([int(line.startswith('TPFAIL')),p,len(B),*B])
    assert [sum(row[0]==t for row in rows) for t in (0,1)]==[5430,126]
    # These auxiliary integer majorants also fit in signed 128 bits at W<=40.
    assert max(2**N*(M+1)**max(N-4,0) for N in range(2,41)
               for M in range(1,41) if N+M-1<=40)<2**127
    data=str(len(rows))+'\n'+'\n'.join(' '.join(map(str,row)) for row in rows)+'\n'
    obj=os.memfd_create('fm151-object');lib=os.memfd_create('fm151-library')
    subprocess.run(['g++','-std=c++17','-O2','-pipe','-fPIC','-fopenmp',
                    '-x','c++','-c','-','-o',f'/proc/self/fd/{obj}'],
                   input=CPP,text=True,pass_fds=(obj,),check=True)
    subprocess.run(['ld','-shared',f'/proc/self/fd/{obj}',
                    '/lib/x86_64-linux-gnu/libstdc++.so.6',
                    '/lib/x86_64-linux-gnu/libgomp.so.1','-lc',
                    '-o',f'/proc/self/fd/{lib}'],
                   pass_fds=(obj,lib),check=True)
    dll=ctypes.CDLL(f'/proc/self/fd/{lib}')
    dll.audit.argtypes=[ctypes.c_char_p];dll.audit.restype=ctypes.c_int
    assert dll.audit(data.encode())==0

def main():
    ap=argparse.ArgumentParser(description='Exact FM-MECH151 checks; memory only.')
    ap.add_argument('--census',help='FM-SEC166 fx5_40.log to audit')
    args=ap.parse_args()
    rng=random.Random(151)
    for _ in range(80):
        C=tuple(rng.choice((-1,1))*rng.randrange(1,9)
                for _ in range(rng.randrange(2,7)))
        p=sum(map(abs,C));d=rng.randrange(max(map(abs,C)),p+1)
        edge(C,d,rng.choice((-1,1)))
    small=0
    for _ in range(100):
        labs=[rng.randrange(1,9) for _ in range(rng.randrange(2,5))]
        signs={n:rng.choice((-1,1)) for n in labs}
        B=tuple(signs[n]*n for n in labs);W=sum(labs)
        for p in range(max(3,max(labs)),W+1):
            if (W-p)%2 or not full_pairfree(B,p): continue
            g=entry(B,p)
            assert g==table(B).get((p,0),0)
            for i in range(len(B)):
                for j in range(i+1,len(B)):
                    if (B[i]-B[j])%2: continue
                    C=B[:i]+B[i+1:j]+B[j+1:]
                    assert g>=entry(C,p);small+=1
    print('boundary checks=80; small comparisons=',small)
    B=(-1,-2,-3,-3,-4,-5,-5,6,8);p=13
    assert check_bounds(B,p)==(8704,149,22540,10934,205,173)
    assert noflip(B,p)[0]==-2
    print('uncovered no-flip example: parent=8704 child=149 budget=-813')
    examples=[
        ((-1,)*42+(2,)*25+(7,25),54,
         2327033494454807709938311397228902169055,237897617235391403865744917550),
        ((-1,)+(2,)*33+(-3,)*21+(-14,),44,
         1284520680193654256810226204082968574465,674771074975898874953071444274164422)]
    for B,p,g,c in examples:
        got=check_bounds(B,p)
        assert got[:2]==(g,c) and g>c
        wrap=(g+2**127)%2**128-2**127
        assert wrap<0
        print('large:',sum(map(abs,B)),len(B),'g=',g,'child=',c,
              'difference=',g-c,'signed128=',wrap)
    coeff=[]
    for j in range(7):
        v=comb(6,j)*76**(6-j)
        if j<=3:v-=64*210*comb(3,j)*82**(3-j)
        if j<=4:v-=4096*comb(4,j)*82**(4-j)
        coeff.append(v)
    assert coeff==[100469760,5908427264,331877376,7422592,82544,456,1]
    assert all(v>0 for v in coeff)
    print('uniform n>=74 polynomial:',coeff)
    if args.census:census(args.census)
    print('ALL CHECKS PASS')

if __name__=='__main__':
    main()
