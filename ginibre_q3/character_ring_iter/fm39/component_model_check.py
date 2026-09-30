# Component model of FM3: F(T) = sum over matchings M of all legs (no chord inside a block) whose
# block-components are each noncrossing diagrams and whose component crossing graph Gamma is bipartite,
# of prod over Gamma-components K=(K1,K2) of (eps^K1 + eps^K2).
import itertools
from collections import defaultdict
def matchings(legs, owner):
    if not legs: yield []; return
    a=legs[0]
    for j in range(1,len(legs)):
        b=legs[j]
        if owner[a]==owner[b]: continue
        for m in matchings(legs[1:j]+legs[j+1:],owner): yield [(a,b)]+m
def cr(e,f):
    a,b=e; c,d=f
    return (a<c<b<d) or (c<a<d<b)
def fusion(labels):
    st={0:1}
    for n in labels:
        nx={}
        for j,c in st.items():
            for k in range(abs(j-n),j+n+1,2): nx[k]=nx.get(k,0)+c
        st=nx
    return st.get(0,0)
def direct(labels,T):
    L=len(labels); tot=0
    for S in range(1<<L):
        a=fusion([labels[i] for i in range(L) if S>>i&1]); b=fusion([labels[i] for i in range(L) if not S>>i&1])
        tot+=(-1)**bin(S&T).count('1')*a*b
    return tot
def model(labels,T):
    owner=[]
    for i,n in enumerate(labels): owner+=[i]*n
    L=len(labels); eps=[-1 if T>>i&1 else 1 for i in range(L)]
    total=0; neg_configs=0; pos_configs=0
    for M in matchings(list(range(len(owner))),owner):
        # block components via chords
        par=list(range(L))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        for a,b in M: par[f(owner[a])]=f(owner[b])
        comp=defaultdict(list)
        for a,b in M: comp[f(owner[a])].append((a,b))
        comps=list(comp.keys())
        # each component noncrossing internally
        ok=all(not cr(e,g) for c in comps for e,g in itertools.combinations(comp[c],2))
        if not ok: continue
        # crossing graph between components
        adj={c:set() for c in comps}
        for c1,c2 in itertools.combinations(comps,2):
            if any(cr(e,g) for e in comp[c1] for g in comp[c2]): adj[c1].add(c2); adj[c2].add(c1)
        # blocks per component
        blocks={c:[i for i in range(L) if f(i)==c] for c in comps}
        # isolated blocks (label 0) none since labels>=1
        color={}; w=1; bip=True
        for c in comps:
            if c in color: continue
            color[c]=0; stack=[c]; side=[[],[]]
            while stack:
                u=stack.pop(); side[color[u]].append(u)
                for v in adj[u]:
                    if v not in color: color[v]=1-color[u]; stack.append(v)
                    elif color[v]==color[u]: bip=False
            e1=1
            for u in side[0]:
                for i in blocks[u]: e1*=eps[i]
            e2=1
            for u in side[1]:
                for i in blocks[u]: e2*=eps[i]
            w*=(e1+e2)
        if not bip: continue
        total+=w
        if w<0: neg_configs+=1
        elif w>0: pos_configs+=1
    return total,neg_configs,pos_configs
bad=0; n=0; stats=[]
for L in range(2,7):
    for lab in itertools.combinations_with_replacement(range(1,4),L):
        if sum(lab)%2 or sum(lab)>12: continue
        for T in range(1<<L):
            if bin(T).count('1')%2: continue
            d=direct(lab,T); m,ng,ps=model(lab,T); n+=1
            if d!=m: bad+=1; print("MISMATCH",lab,T,d,m)
            stats.append((ng,ps,lab,T))
print("checks",n,"mismatches",bad)
print("profiles with negative configurations:",sum(1 for s in stats if s[0]>0))
print("max neg configs:",max(stats)[:4])
