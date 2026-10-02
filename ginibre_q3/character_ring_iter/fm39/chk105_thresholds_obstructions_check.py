from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import comb, factorial

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def mu(ns):
    row = {0: 1}
    for n in ns:
        nxt = {}
        for j, c in row.items():
            for k in cg(j, n):
                nxt[k] = nxt.get(k, 0) + c
        row = nxt
    return row.get(0, 0)

def fusion_row(ns):
    d = {0: 1}
    for n in ns:
        q = {}
        for a, c in d.items():
            for b in cg(a, n):
                q[b] = q.get(b, 0) + c
        d = q
    return d

# Theorem 1: short exhaustive checks.
t1 = 0
for r in range(1, 6):
    for k in range(3, 7):
        for ns in combinations_with_replacement(range(1, 13), k):
            if sum(n >= 2*r for n in ns) >= 3:
                f = fusion_row(ns)
                assert f.get(2*r, 0) >= (r+1) * f.get(0, 0)
                t1 += 1
for r in range(1, 6):
    for a in range(2*r, 2*r+7):
        for b in range(a, 2*r+7):
            for c in range(b, 2*r+7):
                for t in range(13):
                    f = fusion_row((a, b, c, t))
                    assert f.get(2*r, 0) >= (r+1) * f.get(0, 0)
                    t1 += 1
f = fusion_row((1, 10, 11))
assert (f.get(0, 0), f.get(4, 0)) == (1, 2)
print("Theorem 1:", t1, "small exact cases PASS; two-large case", (1, 2))

# Exact parity-class channel profiles and cut-count bound.
def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

@lru_cache(None)
def floor_profile(u, s, v, negative=False):
    ns = [u+(u % 2)]*(s-v) + [u+(1-u % 2)]*v
    if s == 2:
        return (1,)*(max(ns)+1)
    if s == 3:
        m = min(ns)
        S = sum(ns)//2
        return tuple(1+min(m, r, S-r) for r in range(S+1))
    ns.sort(reverse=True)
    hi, third = ns[0], ns[2]
    if s == 4 and negative and min(ns) >= 2:
        hi = max(hi, min(ns)+2)
        base = (Q(1), Q(8,3), Q(10,3), Q(3))
        return tuple(max(r+1 if 2*r <= third else 1,
                         base[r] if r < 4 else Q(0))
                     for r in range(max(hi, 3)+1))
    return tuple(r+1 if 2*r <= third else 1 for r in range(hi+1))

@lru_cache(None)
def cap(u, s, v, t, w, negative=False):
    return sum(a*b for a, b in zip(floor_profile(u,s,v,negative),
                                   floor_profile(u,t,w,negative)))

@lru_cache(None)
def negative_cuts(L, o, q, j, s, v):
    return sum(choose(j,i)*choose(o-j,v-i)*
               choose(q-j,h)*choose(L-o-q+j,s-v-h)
               for i in range(v+1) for h in range(s-v+1)
               if (i+h) % 2)

def threshold_margin(u, L):
    worst = Q(0)
    for o in range(0, L+1, 2):
        childbounds = []
        for rp in (0, 1):
            if (o if rp else L-o) < 2:
                continue
            oc = o-2*rp
            N = L-2
            P = Q(1)
            for s in range(2, N//2+1):
                for v in range(0, s+1, 2):
                    if v > oc or s-v > N-oc:
                        continue
                    P += Q(choose(oc,v)*choose(N-oc,s-v),
                           cap(u,s,v,N-s,oc-v)*(2 if 2*s == N else 1))
            childbounds.append(P/(u+(rp-u) % 2+1))
        if not childbounds:
            continue
        alpha = max(childbounds)
        for q in range(0, L+1, 2):
            for j in range(max(0,o+q-L), min(o,q)+1):
                beta = D = Q(0)
                for s in range(3, L//2+1):
                    for v in range(0, s+1, 2):
                        if v > o or s-v > L-o:
                            continue
                        e = negative_cuts(L,o,q,j,s,v)
                        if not e:
                            continue
                        z = Q(e, cap(u,s,v,L-s,o-v,True) *
                              (2 if 2*s == L else 1))
                        beta += z
                        D += z*s*(L-s)
                worst = max(worst, alpha+beta+D/4096)
    return 1-worst

thresholds = {7:3, 8:7, 9:11, 10:15, 11:21}
expected = {
    7:Q(106131,1713920), 8:Q(135652459,1597736448),
    9:Q(640107145,3846941696), 10:Q(10856088501515,89337731741696),
    11:Q(27696112712089,139929209853440)
}
for L, u in thresholds.items():
    margin = threshold_margin(u, L)
    assert margin == expected[L] and margin > Q(1,20)
    print("threshold margin", L, u, margin)

# Actual signed removal descent on bounded label windows.
def fwht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for lo in range(0, len(v), 2*h):
            for i in range(lo, lo+h):
                a, b = v[i], v[i+h]
                v[i], v[i+h] = a+b, a-b
        h *= 2
    return v

def g_row(B, p):
    n = len(B)
    return fwht([mu(tuple(B[i] for i in range(n) if s >> i & 1)) *
                 mu((p,) + tuple(B[i] for i in range(n)
                                  if not (s >> i & 1)))
                 for s in range(1 << n)])

def sign_classes(xs):
    out = []
    for i, x in enumerate(xs):
        if i == 0 or xs[i-1] != x:
            out.append(0)
        out[-1] |= 1 << i
    return out

for L, u in ((7,3), (8,7), (9,11)):
    profiles = signings = checks = positive_checks = 0
    min_drop = min_slack = None
    for full in combinations_with_replacement(range(u,u+5), L):
        profiles += 1
        B, p = full[:-1], full[-1]
        M = mu(full)
        parent = g_row(B, p)
        cls = sign_classes(B)
        pairs = [(i,j) for i,j in combinations(range(L-1),2)
                 if (B[i]-B[j]) % 2 == 0]
        children = {}
        for i,j in pairs:
            cb = tuple(x for k,x in enumerate(B) if k not in (i,j))
            children[i,j] = (cb, g_row(cb,p))
        pclass = next((c for c in cls if B[c.bit_length()-1] == p), 0)
        for choice in range(1 << len(cls)):
            sg = 0
            for j,c in enumerate(cls):
                if choice >> j & 1:
                    sg |= c
            # The distinguished sign is the product of background signs;
            # enforce pair-freeness if its label occurs in B.
            if pclass and ((sg.bit_count() & 1) != bool(sg & pclass)):
                continue
            signings += 1
            gp = parent[sg]
            for (i,j),(cb,cv) in children.items():
                csg = z = 0
                for k in range(L-1):
                    if k in (i,j):
                        continue
                    if sg >> k & 1:
                        csg |= 1 << z
                    z += 1
                drop = gp-cv[csg]
                slack = 20*drop-M
                assert slack >= 0
                checks += 1
                if M > 0:
                    positive_checks += 1
                    min_drop = drop if min_drop is None else min(min_drop,drop)
                    min_slack = slack if min_slack is None else min(min_slack,slack)
    print("window",L,u,"..",u+4,"profiles/signings/removals",
          profiles,signings,checks,"positive-M checks",positive_checks,
          "minimum drop/slack",min_drop,min_slack)

# Exact polynomial certificates for the concrete FM-MECH162 family.
def char_factor(n, eps):
    out = {}
    for h in range(n//2+1):
        d = n-2*h
        base = Q((-1)**h*choose(n-h,h), 2**d)
        for j in range(d+1):
            z = base*choose(d,j)*(1+eps*(-1)**j)
            if z:
                out[d-j,j] = out.get((d-j,j), Q(0))+z
    return out

def poly_mul(P, R):
    out = {}
    for (i,j), a in P.items():
        for (k,l), b in R.items():
            out[i+k,j+l] = out.get((i+k,j+l), Q(0))+a*b
    return {k:v for k,v in out.items() if v}

def word_poly(factors):
    P = {(0,0):Q(1)}
    for n,e in factors:
        P = poly_mul(P,char_factor(n,e))
    assert all(i%2 == 0 and j%2 == 0 for i,j in P)
    return {(i//2,j//2):v for (i,j),v in P.items() if v}

def padd(P, R):
    z = [Q(0)]*max(len(P),len(R))
    for i,a in enumerate(P): z[i] += a
    for i,a in enumerate(R): z[i] += a
    while len(z)>1 and z[-1] == 0: z.pop()
    return z

def pscale(P, a):
    return [a*x for x in P]

def linear_factor(P, c):
    z = [Q(0)]*(len(P)+1)
    for i,a in enumerate(P):
        z[i] += c*a
        z[i+1] += a
    return z

def rising(start, n):
    P = [Q(1)]
    for k in range(n):
        P = linear_factor(P,start+k)
    return P

def numerator(P, D):
    out = [Q(0)]
    for (i,j),a in P.items():
        term = [a*16**i*Q(factorial(2*j)*factorial(2*j+1),
                           factorial(j)**2)]
        for start,n in ((Q(1,2),i),(Q(3,2),i),
                        (Q(2+i+j),D-i-j),(Q(3+i+j),D-i-j)):
            q = rising(start,n)
            z = [Q(0)]*(len(term)+len(q)-1)
            for x,a0 in enumerate(term):
                for y,b0 in enumerate(q):
                    z[x+y] += a0*b0
            term = z
        out = padd(out,term)
    return out

def shift_poly(P, s):
    q = [Q(0)]*len(P)
    for d,a in enumerate(P):
        for k in range(d+1):
            q[k] += a*choose(d,k)*s**(d-k)
    return q

core = [(2,-1),(4,-1),(6,-1),(8,-1),(10,1)]
P = word_poly(core)
D = max(i+j for i,j in P)
base = numerator(P,D)
child = word_poly([(2,-1),(4,-1),(10,1)])
certs = [(base,1024),(padd(numerator(child,D),pscale(base,-1)),1024)]
for i,j in combinations(range(5),2):
    z = [list(x) for x in core]
    z[i][1] *= -1
    z[j][1] *= -1
    certs.append((padd(numerator(word_poly([tuple(x) for x in z]),D),
                       pscale(base,-1)),1024))
plus_pair = {}
for (i,j),a in P.items():
    plus_pair[i+1,j] = plus_pair.get((i+1,j),Q(0))+a
    plus_pair[i,j+1] = plus_pair.get((i,j+1),Q(0))-a
certs.append((numerator(plus_pair,D+1),1023))
assert len(certs) == 13
assert sum(len(shift_poly(P,s)) for P,s in certs) == 367
assert all(min(shift_poly(P,s)) > 0 for P,s in certs)
print("Theorem 3: 13 shifted polynomial certificates, 367 positive coefficients PASS")

# Proposition 4: actual signed table and the artificial full-coordinate change.
def g_value(B,p):
    z = 0
    for s in range(1 << len(B)):
        a = tuple(sorted(abs(B[i]) for i in range(len(B)) if s >> i & 1))
        b = (p,) + tuple(sorted(abs(B[i]) for i in range(len(B))
                                if not (s >> i & 1)))
        z += (-1 if s.bit_count() & 1 else 1)*mu(a)*mu(b)
    return z

def walsh(sel, minus, m):
    z = 0
    s = sel
    while True:
        z += (-1 if (s & minus).bit_count() & 1 else 1)*m[s]*m[sel^s]
        if not s:
            break
        s = (s-1) & sel
    return z

B = (-1,-3,-3,-4,-5,-5,-6)
p = 7
L = 8
full = (1 << L)-1
labels = tuple(abs(x) for x in B)+(p,)
m = [mu(tuple(sorted(labels[i] for i in range(L) if s >> i & 1)))
     for s in range(1 << L)]
neg = sum(1 << i for i,x in enumerate(B) if x < 0)
weights = [(-1 if (s & neg).bit_count() & 1 else 1)*m[s]*m[full^s]
           for s in range(1 << 7)]
pairs = [(i,j) for i,j in combinations(range(7),2)
         if (B[i]-B[j]) % 2 == 0]
I,J = max(pairs,key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                               max(abs(B[ij[0]]),abs(B[ij[1]]))))
assert (I,J) == (3,6)
blocks = [0,0,0,0]
for s,z in enumerate(weights):
    blocks[((s >> I)&1)+2*((s >> J)&1)] += z
T,X,Y,Z = blocks[0],blocks[2],blocks[1],blocks[3]
child_B = tuple(B[k] for k in range(7) if k not in (I,J))
h = g_value(child_B,p)
assert (m[full],T,X,Y,Z,sum(weights),h) == (624,741,-98,-164,71,550,31)

cuts = []
for i,j in combinations(range(L),2):
    cuts.append(sum(weights[s] for s in range(1 << 7)
        if ((s >> i)&1 if i < 7 else 0) !=
           ((s >> j)&1 if j < 7 else 0)))
assert len(cuts) == 28 and max(cuts) == -3
proper_products = [m[s]*m[full^s] for s in range(1,full)]
assert max(proper_products) == 46

proper_count = 0
for sel in range(full):
    q = sel
    while True:
        if q.bit_count() % 2 == 0:
            assert walsh(sel,q,m) >= 0
            proper_count += 1
        if not q:
            break
        q = (q-1) & sel
assert proper_count == 3153
assert min(walsh(full,q,m) for q in range(1 << L)
           if q.bit_count() % 2 == 0) == 1060

modified = m.copy()
modified[full] = 104
assert modified[full] >= max(proper_products)
assert min(walsh(full,q,modified) for q in range(1 << L)
           if q.bit_count() % 2 == 0) == 20
assert walsh(full,full,modified)//2 == 30 < h
assert T-520-h-abs(X)-abs(Y)-abs(Z) == -143
print("Proposition 4: 28 cuts max -3; 3153 proper patterns PASS;",
      "split max 46; full minima 1060/20; parent 30 < child 31; S=-143")
print("FM-CHK105 standalone checks PASS")
