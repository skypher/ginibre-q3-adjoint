"""FM-SEC109 (luna_max_jupiter): Hall/max-flow with intermediate admissible targets (Gamma*); first failure (1,1,1,1,2,2), T={1,2}."""
from functools import lru_cache
from itertools import combinations, product
from collections import defaultdict, deque

def cross(e,f):
    a,b=sorted(e); c,d=sorted(f)
    return a<c<b<d or c<a<d<b

def smooth(M,e,f):
    a,b,c,d=sorted(tuple(e)+tuple(f))
    B=set(M)-{tuple(sorted(e)),tuple(sorted(f))}
    return (tuple(sorted(B|{(a,b),(c,d)})),
            tuple(sorted(B|{(a,d),(b,c)})))

def owners(ns):
    return tuple(i for i,n in enumerate(ns) for _ in range(n))

@lru_cache(None)
def matchings(ns):
    own=owners(ns)
    @lru_cache(None)
    def rec(rem):
        if not rem:
            return ((),)
        a=rem[0]
        out=[]
        for j in range(1,len(rem)):
            b=rem[j]
            if own[a]==own[b]:
                continue
            for tail in rec(rem[1:j]+rem[j+1:]):
                out.append(tuple(sorted(((a,b),)+tail)))
        return tuple(out)
    return rec(tuple(range(sum(ns))))

@lru_cache(None)
def factors(ns,M):
    own=owners(ns)
    parent=list(range(len(ns)))

    def find(x):
        if parent[x]!=x:
            parent[x]=find(parent[x])
        return parent[x]

    for a,b in M:
        if own[a]==own[b]:
            return None
        x,y=find(own[a]),find(own[b])
        if x!=y:
            parent[x]=y

    groups=defaultdict(list)
    for e in M:
        groups[find(own[e[0]])].append(e)
    E=list(groups.values())

    if any(cross(e,f) for G in E for e,f in combinations(G,2)):
        return None

    masks=[]
    for G in E:
        mask=0
        for e in G:
            for v in e:
                mask |= 1 << own[v]
        masks.append(mask)

    adj=[set() for _ in E]
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            if any(cross(e,f) for e in E[i] for f in E[j]):
                adj[i].add(j)
                adj[j].add(i)

    color={}
    out=[]
    for s in range(len(E)):
        if s in color:
            continue
        color[s]=0
        todo=[s]
        side=[0,0]
        while todo:
            u=todo.pop()
            side[color[u]] |= masks[u]
            for v in adj[u]:
                if v not in color:
                    color[v]=1-color[u]
                    todo.append(v)
                elif color[v]==color[u]:
                    return None
        out.append(tuple(side))
    return tuple(out)

def weight(F,T):
    if F is None:
        return None
    ans=1
    for x,y in F:
        ex=-1 if (x&T).bit_count()%2 else 1
        ey=-1 if (y&T).bit_count()%2 else 1
        ans *= ex+ey
    return ans

@lru_cache(None)
def reachable(ns,M,T=None,stop_zero=False):
    own=owners(ns)
    seen={M}
    stack=[M]
    while stack:
        X=stack.pop()
        if stop_zero and X!=M:
            F=factors(ns,X)
            if F is not None and weight(F,T)==0:
                continue
        for e,f in combinations(X,2):
            if not cross(e,f):
                continue
            for Y in smooth(X,e,f):
                if any(own[a]==own[b] for a,b in Y):
                    continue
                if Y not in seen:
                    seen.add(Y)
                    stack.append(Y)
    return frozenset(seen-{M})

def maxflow(demand,capacity,edges):
    n,m=len(demand),len(capacity)
    source=n+m
    sink=source+1
    graph=[[] for _ in range(sink+1)]

    def add(u,v,c):
        graph[u].append([v,c,len(graph[v])])
        graph[v].append([u,0,len(graph[u])-1])

    inf=sum(demand)
    for i,x in enumerate(demand):
        add(source,i,x)
    for j,x in enumerate(capacity):
        add(n+j,sink,x)
    for i,js in enumerate(edges):
        for j in js:
            add(i,n+j,inf)

    flow=0
    while True:
        level=[-1]*len(graph)
        level[source]=0
        q=deque([source])
        while q:
            u=q.popleft()
            for v,c,_ in graph[u]:
                if c and level[v]<0:
                    level[v]=level[u]+1
                    q.append(v)
        if level[sink]<0:
            return flow

        it=[0]*len(graph)
        def dfs(u,f):
            if u==sink:
                return f
            while it[u]<len(graph[u]):
                edge=graph[u][it[u]]
                v,c,rev=edge
                if c and level[v]==level[u]+1:
                    sent=dfs(v,min(f,c))
                    if sent:
                        edge[1]-=sent
                        graph[v][rev][1]+=sent
                        return sent
                it[u]+=1
            return 0

        while True:
            sent=dfs(source,inf)
            if not sent:
                break
            flow+=sent

def analyze(ns,T,stop_zero=False):
    Ms=matchings(ns)
    FF={M:factors(ns,M) for M in Ms}
    neg=[]
    pos=[]
    demand=[]
    capacity=[]
    for M in Ms:
        w=weight(FF[M],T)
        if w is not None and w<0:
            neg.append(M)
            demand.append(-w)
        elif w is not None and w>0:
            pos.append(M)
            capacity.append(w)

    ix={M:j for j,M in enumerate(pos)}
    edges=[
        {ix[N] for N in reachable(ns,M,T,stop_zero) if N in ix}
        for M in neg
    ]
    flow=maxflow(demand,capacity,edges) if neg else 0
    return {
        'Ms':Ms, 'FF':FF, 'neg':neg, 'pos':pos,
        'demand':demand, 'capacity':capacity, 'edges':edges,
        'flow':flow,
        'F':sum(weight(FF[M],T) or 0 for M in Ms)
    }

profiles=[]
for L in range(1,8):
    bound=4 if L<=6 else 3
    for ns in product(range(1,bound+1),repeat=L):
        S=sum(ns)
        if S<=8 and S%2==0 and max(ns)<=S//2:
            profiles.append(ns)
profiles.sort(key=lambda ns:(sum(ns),len(ns),ns))

checked=0
for ns in profiles:
    for T in range(1<<len(ns)):
        if T.bit_count()%2:
            continue
        checked+=1
        z=analyze(ns,T)
        if z['flow']<sum(z['demand']):
            print('FIRST',ns,
                  'T',tuple(i for i in range(len(ns)) if T>>i&1),
                  'flow/demand',z['flow'],sum(z['demand']),
                  'positive capacity',sum(z['capacity']),
                  'F',z['F'],'pairs checked',checked)
            break
    else:
        continue
    break

ns=(1,1,1,1,2,2)
T=(1<<1)|(1<<2)
z=analyze(ns,T)
z0=analyze(ns,T,stop_zero=True)
print('FLOW_FULL_AND_ZERO_BLOCKED',
      z['flow'],z0['flow'],'DEMAND',sum(z['demand']))

for i,M in enumerate(z['neg']):
    print('SOURCE',i,M,'weight',-z['demand'][i],
          'targets',[z['pos'][j] for j in sorted(z['edges'][i])])

hall=[]
for mask in range(1,1<<len(z['neg'])):
    src=[i for i in range(len(z['neg'])) if mask>>i&1]
    Nset=set().union(*(z['edges'][i] for i in src))
    d=sum(z['demand'][i] for i in src)
    c=sum(z['capacity'][j] for j in Nset)
    hall.append((c-d,tuple(src),d,c,tuple(sorted(Nset))))
print('MIN_HALL_SLACK',min(hall))
print('POSITIVE_TARGETS',
      [(N,z['capacity'][j]) for j,N in enumerate(z['pos'])])

for ns,T in [
    ((1,1,1,1,2),5),
    ((2,2,2,2,2,2),48),
    ((1,1,1,1,1,1,1,1,6),15)
]:
    y=analyze(ns,T)
    print('BOUNDARY',ns,T,
          'flow/demand',y['flow'],sum(y['demand']),
          'capacity',sum(y['capacity']),'F',y['F'])
