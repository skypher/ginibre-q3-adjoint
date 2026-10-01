# Adversarial local search for small g_p / m_p(all) over PAIR-FREE unsaturated lists with k >= 2 non-distinguished cores.
import random, sys, time
from collections import defaultdict
from multiprocessing import Pool
def gv(lab):
    st={(0,0):1}
    for x in lab:
        n=abs(x); e=1 if x>0 else -1; nx=defaultdict(int)
        for (s,t),v in st.items():
            for s2 in range(abs(s-n),s+n+1,2): nx[s2,t]+=v
            for t2 in range(abs(t-n),t+n+1,2): nx[s,t2]+=e*v
        st=nx
    return st
def mv(lab):
    st={0:1}
    for x in lab:
        n=abs(x); nx=defaultdict(int)
        for s,v in st.items():
            for s2 in range(abs(s-n),s+n+1,2): nx[s2]+=v
        st=nx
    return st
def score(bg):
    # best (smallest) ratio over admissible p for this background
    if any((v in bg) and (-v in bg) for v in set(bg)): return None
    W=sum(abs(v) for v in bg); mx=max(abs(v) for v in bg)
    if sum(1 for v in bg if abs(v)>=3) < 2: return None
    g=gv(bg); m=mv(bg); nminus=sum(1 for v in bg if v<0); psign=1 if nminus%2==0 else -1
    best=None
    for p in range(max(mx,3), W+1):
        if (W-p)%2: continue
        d=(W-p)//2
        if d<8 or mx>d: continue
        if (-psign*p) in bg: continue           # distinguished must keep the list pair-free
        G=g.get((p,0),0); M=m.get(p,0)
        if G<0: return (-1e9, p, G)
        if M>0:
            r=G/M
            if best is None or r<best[0]: best=(r,p,G)
    return best
def mutate(bg, rng):
    b=list(bg); op=rng.random()
    vals=sorted(set(abs(v) for v in b))
    if op<0.25 and len(b)>3: b.pop(rng.randrange(len(b)))
    elif op<0.5:
        v=rng.choice(b); b.append(v)                          # duplicate an existing factor (same sign)
    elif op<0.7:
        n=rng.randint(1,20); sg=rng.choice([1,-1])
        if (-sg*n) in b: sg=-sg
        b.append(sg*n)
    elif op<0.85:
        n=rng.choice(vals); b=[(-v if abs(v)==n else v) for v in b]   # flip the sign of a whole label value
    else:
        i=rng.randrange(len(b)); n=abs(b[i]); sg=1 if b[i]>0 else -1; n2=max(1,n+rng.choice([-1,1]))
        b[i]=sg*n2
        if (-sg*n2) in b: return bg
    return b
def run(seed):
    rng=random.Random(seed); best_global=None; t0=time.time()
    while time.time()-t0 < float(sys.argv[1]):
        k=rng.randint(2,4); bg=[rng.choice([1,-1])*rng.randint(3,12) for _ in range(k)]
        s1=rng.choice([1,-1]); bg+= [s1]*rng.randint(4,24); bg+=[rng.choice([2,-2])]*rng.randint(0,6)
        # make pair-free: unify signs per value
        sgn={}
        bg=[(sgn.setdefault(abs(v), 1 if v>0 else -1))*abs(v) for v in bg]
        cur=score(bg)
        if cur is None: continue
        for it in range(150):
            nb=mutate(bg, rng)
            if sum(abs(v) for v in nb) > WMAX: continue
            sc=score(nb)
            if sc is None: continue
            if sc[0] <= cur[0]: bg, cur = nb, sc
        if best_global is None or cur[0] < best_global[0]: best_global=(cur[0], cur[1], cur[2], sorted(bg))
    return best_global
WMAX=int(sys.argv[2]) if len(sys.argv)>2 else 90
if __name__=="__main__":
    with Pool(32) as pool:
        res=pool.map(run, range(32))
    res=[r for r in res if r]
    res.sort(key=lambda r: r[0])
    for r in res[:8]: print(f"ratio {r[0]:.4f} p {r[1]} g {r[2]} background {r[3]}")
    print("NEGATIVE FOUND" if res and res[0][0] < 0 else "no negative")
