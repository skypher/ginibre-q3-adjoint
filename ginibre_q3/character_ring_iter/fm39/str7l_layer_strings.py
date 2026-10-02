from collections import defaultdict
from datetime import datetime, timezone

def stamp(s):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), s, flush=True)
def cg(a,b): return range(abs(a-b), a+b+1, 2)
def mul(P,Q):
    out=defaultdict(int)
    for (a,b),u in P.items():
        for (c,d),v in Q.items():
            for x in cg(a,c):
                for y in cg(b,d): out[x,y]+=u*v
    return {k:v for k,v in out.items() if v}
def dmul(P):
    out=defaultdict(int)
    for (i,j),v in P.items():
        out[i+1,j]+=v
        if i: out[i-1,j]+=v
        out[i,j+1]-=v
        if j: out[i,j-1]-=v
    return {k:v for k,v in out.items() if v}
def quotient(word):
    P={(0,0):1}
    for z in word:
        n=abs(z)
        F=({(j,n-1-j):1 for j in range(n)} if z<0
           else {(n,0):1,(0,n):1})
        P=mul(P,F)
    return P
def sp4_vector(G):
    dG=dmul(G)
    assert all(dG.get((j,i),0)==-v for (i,j),v in dG.items())
    assert all(v==0 for (i,j),v in dG.items() if i==j)
    v={(i-1,j):c for (i,j),c in dG.items() if i>j and c}
    rebuild=defaultdict(int)
    for (a,b),c in v.items():
        rebuild[a+1,b]+=c
        rebuild[b,a+1]-=c
    assert {k:c for k,c in rebuild.items() if c}==dG
    return v
def signed_table(word):
    P={(0,0):1}
    for z in word:
        n=abs(z); eps=1 if z>0 else -1
        out=defaultdict(int)
        for (a,b),v in P.items():
            for c in cg(a,n): out[c,b]+=v
            for c in cg(b,n): out[a,c]+=eps*v
        P={k:v for k,v in out.items() if v}
    return P
def interior_cut(C):
    A=[]; B=[]; wa=wb=0
    for z in sorted(C,key=lambda z:(-abs(z),z)):
        if wa<=wb: A.append(z); wa+=abs(z)
        else: B.append(z); wb+=abs(z)
    if wa>wb: A,B=B,A
    return A,B
def l_image(lam,T):
    a,b=lam; ell=a-b
    if T<ell: return None
    j=min(b,(T-ell)//2)
    return ell+j,j
def rows(V):
    qs=sorted({a+b for a,b in V})
    return {q:tuple(V.get((q-b,b),0) for b in range(q//2+1)) for q in qs}
def run_case(L,ij,name,allocation):
    i,j=ij; C=[z for k,z in enumerate(L) if k not in (i,j)]
    A,B=interior_cut(C); Bp=B+[L[i],L[j]]
    ma=sum(z<0 for z in A); mb=sum(z<0 for z in Bp); ra,rb=allocation
    assert ma>=ra and mb>=rb and (ma-ra)%2==0 and (mb-rb)%2==0
    qa,qb=quotient(A),quotient(Bp); ga,gb=qa,qb
    for _ in range(ma-ra): ga=dmul(ga)
    for _ in range(mb-rb): gb=dmul(gb)
    va,vb=sp4_vector(ga),sp4_vector(gb)
    fa,fb=signed_table(A),signed_table(Bp); layers=defaultdict(int)
    for (r,s),x in fa.items(): layers[r+s]+=x*fb.get((r,s),0)
    net=defaultdict(int); pos=defaultdict(int); neg=defaultdict(int)
    for lam,x in va.items():
        p=x*vb.get(lam,0); q=sum(lam); net[q]+=p
        if p>=0: pos[q]+=p
        else: neg[q]-=p
    def gamma(T):
        if allocation==(1,1):
            return sum(x*vb.get(lam,0) for lam,x in va.items() if sum(lam)<=T-1)
        assert allocation==(2,0)
        return sum(y*va.get(l_image(lam,T),0) for lam,y in vb.items()
                   if l_image(lam,T) is not None)
    for T in range(-1,sum(map(abs,L))+3):
        assert sum(x for t,x in layers.items() if t<=T)==2*gamma(T),(name,T)
    for t in range(sum(map(abs,L))+3):
        assert layers.get(t,0)==2*(gamma(t)-gamma(t-1)),(name,t)
    negA=sorted((lam,x) for lam,x in va.items() if x<0)
    stamp(f'{name}: A={A}; Bprime={Bp}; minus=({ma},{mb}); allocation={allocation}; vector_sizes=({len(va)},{len(vb)}); negative_counts=({len(negA)},{sum(x<0 for x in vb.values())}); minima=({min(va.values())},{min(vb.values())})')
    if name!='minus4,plus6':
        stamp(f'{name}: full GA rows (q; beta=0..floor(q/2))={rows(va)}')
        stamp(f'{name}: full GB rows (q; beta=0..floor(q/2))={rows(vb)}')
    stamp(f'{name}: GA negative entries={negA}')
    stamp(f'{name}: q rows (q, positive products, negative products, net)={[(q,pos[q],neg[q],net[q]) for q in sorted(net)]}')
    nonzero=sorted((t,x) for t,x in layers.items() if x)
    stamp(f'{name}: layer_min={min(nonzero,key=lambda z:z[1])}; negative_layers={[(t,x) for t,x in nonzero if x<0]}; Phi={sum(layers.values())}')
    return va,vb,dict(net),dict(pos),dict(neg),dict(layers)

L=[1,1,-2]+[3]*6+[-4]*3+[5]*4+[6]
stamp('FM-STR7l exact W1 verifier start; W=60, 17 factors')
mm=run_case(L,(9,10),'minus4,minus4',(2,0))
m5=run_case(L,(9,12),'minus4,plus5',(2,0))
m6=run_case(L,(9,16),'minus4,plus6',(2,0))
top=run_case(L,(12,13),'TopPair plus5,plus5',(1,1))
assert mm[2]=={0:1332240,2:4457760,4:9613999,6:12589428,8:10927633,10:9508230,12:7189912,14:2999020,16:1170115,18:465268,20:19539,22:-10330,24:7961}
assert mm[5][22]==432908 and all(x>=0 for x in mm[5].values())
assert m5[2]==m6[2]=={1:3346606,3:8018080,5:10613013,7:12580221,9:11550374,11:7088412,13:4057562,15:2280227,17:587988,19:80029,21:57435,23:10828}
assert all(x>=0 for x in m5[5].values()) and m5[0]==m6[0] and m5[1]==m6[1]
assert m5[0].get((11,10))==m5[0].get((13,10))==-3
assert m5[1].get((11,10))==-2457 and m5[1].get((13,10))==-1442
assert top[2]=={0:1969016,2:6240746,4:13772427,6:15776316,8:13978391,10:7657964,12:3782102,14:-461364,16:-1202647,18:-711065,20:-409610,22:-114696,24:-6805}
assert {t:x for t,x in top[5].items() if x<0}=={15:-922728,17:-2405294,19:-1422130,21:-819220,23:-229392,25:-13610}
assert top[0].get((16,16),0)==0 and top[1][(16,16)]==-12
assert sum(mm[5].values())==sum(m5[5].values())==sum(m6[5].values())==sum(top[5].values())==120541550
stamp('PASS: full Sp(4) vectors, size-slice products, exact L_T prefixes, and raw layers')
