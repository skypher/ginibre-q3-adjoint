from collections import Counter, defaultdict
from itertools import product
from random import Random
from datetime import datetime, timezone
from fractions import Fraction
import re
import time

import gzip as _gz, pathlib as _pl
LOG = _pl.Path(__file__).resolve().parent / "sec166_census_w52_noflip.log.gz"
DIRS = [(1,1), (1,0), (0,1), (2,1), (1,2), (2,-1), (-1,2)]

def stamp():
    return datetime.now(timezone.utc).strftime("%H:%M:%S UTC")

def fusion(i, j):
    return range(abs(i-j), i+j+1, 2)

def character(xs):
    T = {(0,0): 1}
    for z in xs:
        n, eps = abs(z), (1 if z > 0 else -1)
        R = defaultdict(int)
        for (i,j), c in T.items():
            for r in fusion(i,n):
                R[r,j] += c
            for s in fusion(j,n):
                R[i,s] += eps*c
        T = {k:v for k,v in R.items() if v}
    return T

def pair_terms(u, v):
    a, b = abs(u), abs(v)
    eu, ev = (1 if u > 0 else -1), (1 if v > 0 else -1)
    pure, mixed = defaultdict(int), defaultdict(int)
    for q in fusion(a,b):
        pure[q,0] += 1
        pure[0,q] += eu*ev
    mixed[a,b] += ev
    mixed[b,a] += eu
    return dict(pure), dict(mixed)

def dot_expand(A, B, G):
    """Height coefficients of sum A[r,s] [B*G][r,s]."""
    rows = {r for r,s in A}
    keys = set(A)
    out = defaultdict(int)
    for (i,j), cb in B.items():
        for (k,l), cg in G.items():
            for r in fusion(i,k):
                if r not in rows:
                    continue
                for s in fusion(j,l):
                    if (r,s) in keys:
                        out[r+s] += A[r,s]*cb*cg
    return dict(out)

def contributions(Afactors, Bfactors, u, v):
    A, B = character(Afactors), character(Bfactors)
    pure, mixed = pair_terms(u,v)
    return dot_expand(A,B,pure), dot_expand(A,B,mixed)

def interior_cut(C):
    A, B, wa, wb = [], [], 0, 0
    for z in sorted(C, key=lambda q:(-abs(q), q)):
        if wa <= wb:
            A.append(z); wa += abs(z)
        else:
            B.append(z); wb += abs(z)
    if wa > wb:
        A, B = B, A
    return A, B

def pair_budget(xs, u, v, A=None, B=None):
    C = list(xs)
    C.remove(u); C.remove(v)
    if A is None:
        A, B = interior_cut(C)
    P, M = contributions(A,B,u,v)
    phi = character(xs).get((0,0),0)
    assert sum(P.values()) + sum(M.values()) == phi
    pref, minimum = 0, 0
    rows = []
    for t in range(max(set(P)|set(M), default=-1)+1):
        p, m = P.get(t,0), M.get(t,0)
        pref += p + min(m,0)
        minimum = min(minimum,pref)
        rows.append((t,p,m,pref))
    return A, B, phi, P, M, rows, minimum

def pair_types(xs):
    c = Counter(xs)
    vals = sorted(c, key=lambda z:(abs(z),z))
    for i,u in enumerate(vals):
        for j in range(i,len(vals)):
            v = vals[j]
            if i == j:
                if c[u] >= 2:
                    yield u,v,c[u]*(c[u]-1)//2
            else:
                yield u,v,c[u]*c[v]

def count_vectors(total, width=5):
    def rec(i, left, pref):
        if i == width-1:
            yield tuple(pref+[left])
            return
        for q in range(left+1):
            yield from rec(i+1,left-q,pref+[q])
    yield from rec(0,total,[])

def profiles_of_length(length, max_label=5):
    for counts in count_vectors(length,max_label):
        labels = [i+1 for i,c in enumerate(counts) if c]
        mults = [counts[n-1] for n in labels]
        for signs in product((-1,1), repeat=len(labels)):
            if sum(m for m,e in zip(mults,signs) if e < 0) % 2:
                continue
            xs = []
            for n,m,e in zip(labels,mults,signs):
                xs.extend([e*n]*m)
            yield tuple(xs)

def prefix_nonnegative(P, M):
    q = 0
    for t in range(max(set(P)|set(M), default=-1)+1):
        q += P.get(t,0) + min(M.get(t,0),0)
        if q < 0:
            return False, (t,q)
    return True, None

def small_all_pair_screen():
    bylen = {}
    profiles = physical = pairtypes = zero_phi = 0
    for L in range(2,11):
        nprof = 0
        for xs in profiles_of_length(L):
            nprof += 1; profiles += 1
            phi = character(xs).get((0,0),0)
            zero_phi += (phi == 0)
            for u,v,mult in pair_types(xs):
                pairtypes += 1; physical += mult
                C = list(xs); C.remove(u); C.remove(v)
                A,B = interior_cut(C)
                P,M = contributions(A,B,u,v)
                assert sum(P.values())+sum(M.values()) == phi
                ok,bad = prefix_nonnegative(P,M)
                if not ok:
                    raise AssertionError(("small pair failure",xs,u,v,bad))
        bylen[L] = nprof
        print(stamp(),"small-pair",L,nprof,"so far",profiles,physical,flush=True)
    print("SMALL_PAIR_RESULT",bylen,profiles,pairtypes,physical,zero_phi)

def all_split_test(A, B):
    ta, tb = character(A), character(B)
    dots = {k:v*tb.get(k,0) for k,v in ta.items() if v*tb.get(k,0)}
    for l1,l2 in DIRS:
        bins = defaultdict(int)
        for (r,s),v in dots.items():
            bins[l1*r+l2*s] += v
        acc = 0
        for T in sorted(bins):
            acc += bins[T]
            if acc < 0:
                return False,(l1,l2,T,acc)
    return True,None

def exhaustive_all_splits():
    nprof = nsplits = checks = 0
    for L in range(2,7):
        for xs in profiles_of_length(L):
            nprof += 1
            for mask in range(1,(1<<L)-1):
                A = [xs[i] for i in range(L) if mask>>i&1]
                B = [xs[i] for i in range(L) if not(mask>>i&1)]
                ok,bad = all_split_test(A,B)
                if not ok:
                    raise AssertionError(("all-split failure",xs,mask,bad))
                nsplits += 1; checks += len(DIRS)
    print("ALL_SPLITS_EXACT",nprof,nsplits,checks,"fail=0")

def random_all_splits():
    rng = Random(1792026)
    checks, maxweight = 0, 0
    for _ in range(600):
        L = rng.randint(4,11)
        labs = [rng.randint(1,18) for _ in range(L)]
        signs = [rng.choice((-1,1)) for _ in range(L)]
        if sum(e<0 for e in signs)%2:
            signs[-1] *= -1
        xs = [e*n for e,n in zip(signs,labs)]
        mask = rng.randrange(1,(1<<L)-1)
        A = [xs[i] for i in range(L) if mask>>i&1]
        B = [xs[i] for i in range(L) if not(mask>>i&1)]
        ok,bad = all_split_test(A,B)
        if not ok:
            raise AssertionError(("random split failure",xs,mask,bad))
        checks += len(DIRS); maxweight=max(maxweight,sum(map(abs,xs)))
    print("ALL_SPLITS_RANDOM",600,checks,"maxweight",maxweight,"fail=0")

def all_cuts_small():
    nprof = pairtypes = cuts = 0
    for L in range(2,7):
        for xs in profiles_of_length(L):
            nprof += 1
            for u,v,_mult in pair_types(xs):
                pairtypes += 1
                C = list(xs); C.remove(u); C.remove(v)
                for mask in range(1<<len(C)):
                    A = [C[i] for i in range(len(C)) if mask>>i&1]
                    B = [C[i] for i in range(len(C)) if not(mask>>i&1)]
                    P,M = contributions(A,B,u,v)
                    ok,bad = prefix_nonnegative(P,M)
                    if not ok:
                        raise AssertionError(("cut failure",xs,u,v,A,B,bad))
                    cuts += 1
    print("ALL_CUTS_SMALL",nprof,pairtypes,cuts,"fail=0")

def top_pair(B):
    candidates = []
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if abs(B[i])%2 == abs(B[j])%2:
                key = (abs(B[i])+abs(B[j]),
                       max(abs(B[i]),abs(B[j])),
                       min(abs(B[i]),abs(B[j])))
                candidates.append((key,B[i],B[j]))
    if not candidates:
        raise ValueError(("no TopPair",B))
    _,u,v = max(candidates)
    return u,v

def noflip_top_pair_replay():
    bands = {
        "<=40": [0,0,0],
        "41-44": [0,0,0],
        "45-48": [0,0,0],
        "49-52": [0,0,0],
    }
    first_failure = None
    with _gz.open(LOG, 'rt') as f:
        for row,line in enumerate(f,1):
            m = re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+)",line)
            if not m:
                raise ValueError(("bad input line",row,line))
            W,p = int(m.group(1)),int(m.group(2))
            B = [int(z) for z in m.group(3).split()]
            phi_log = int(m.group(4))
            xs = B+[p]
            phi = character(xs).get((0,0),0)
            assert phi == phi_log
            u,v = top_pair(B)
            A,C = list(xs),None
            A.remove(u); A.remove(v)
            # interior_cut receives the full list with the selected pair removed.
            left,right = interior_cut(A)
            P,M = contributions(left,right,u,v)
            assert sum(P.values())+sum(M.values()) == phi
            ok,bad = prefix_nonnegative(P,M)
            band = "<=40" if W<=40 else ("41-44" if W<=44 else ("45-48" if W<=48 else "49-52"))
            bands[band][0] += 1
            bands[band][1] += 1
            if not ok and first_failure is None:
                first_failure = (row,W,p,B,u,v,bad,phi)
            if row%5000==0:
                print(stamp(),"no-flip rows",row,"failure",first_failure,flush=True)
    print("NOFLIP_TOPPAIR_RESULT",bands,"failure",first_failure)

def named_split_sample():
    rng = Random(1792026)
    fams = {
      "even12":[-40,42,-44,46,-48,50,-52,54,-56,58,-60,62],
      "F1":[1,1,-2]+[3]*4+[-4]*3+[5]*4+[8],
    }
    for name,xs in fams.items():
        for rep in range(3):
            ids = set(rng.sample(range(len(xs)),rng.randint(1,len(xs)-1)))
            A = [x for i,x in enumerate(xs) if i in ids]
            B = [x for i,x in enumerate(xs) if i not in ids]
            ok,bad = all_split_test(A,B)
            print("NAMED_SPLIT",name,rep+1,len(A),len(B),ok,bad)
            assert ok

if __name__ == "__main__":
    # Each section prints progress; the log replay uses the supplied W<=52 census.
    exhaustive_all_splits()
    random_all_splits()
    all_cuts_small()
    small_all_pair_screen()
    named_split_sample()
    noflip_top_pair_replay()
