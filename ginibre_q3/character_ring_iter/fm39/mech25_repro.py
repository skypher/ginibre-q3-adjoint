
import argparse
from math import comb
from functools import lru_cache

ap = argparse.ArgumentParser()
ap.add_argument("--center-max", type=int, default=100)
ap.add_argument("--metric-max", type=int, default=90)
args = ap.parse_args()

def row(a, e):
    N, d = a + e, a - e
    c, prev = [1], 0
    for k in range(N):
        num = d*c[-1] - (N-k+1)*prev
        assert num % (k+1) == 0
        prev, nxt = c[-1], num // (k+1)
        c.append(nxt)
    return c

def half(v):
    return 0 if v[1] > 0 or (v[1] == 0 and v[0] > 0) else 1

def lt(v, w):
    if half(v) != half(w):
        return half(v) < half(w)
    return v[0]*w[1] - v[1]*w[0] > 0

def long_sweep(c, x, y):
    N = len(c)-1
    at = lambda k: c[k] if 0 <= k <= N else 0
    g = [(at(k), at(k-1)) for k in range(x, y+2)]
    turns = sum(lt(v, u) for u, v in zip(g, g[1:]))
    return turns >= 2 or (turns == 1 and lt(g[0], g[-1]))

def domain(c, lo, hi):
    for p in range(lo+1, hi):
        if c[p] == 0:
            continue
        sig = 1 if c[p] > 0 else -1
        if sig*c[p-1] > 0:
            continue
        q = p
        while q+1 <= hi and sig*c[q+1] > 0:
            q += 1
        if q < hi:
            return p, q
    return None

def metric(N, d, x, C):
    h, n = 2*x+C-N, C+1
    L = max(1, n-h)
    eta = n*n if 2*L >= n else 4*L*(n-L)
    V = ((N+2)**2-d*d)//4
    R, r = C+h, abs(C-h)
    K = (N+2-C)**2-h*h
    A = eta*(4*V+R*r)-K
    ok = (R*R < 4*V and A >= 0
          and A*A >= 4*eta*eta*(R+r)**2*V)
    return ok, L, eta, V, R, r, K

branches = dict(nonpositive=0, wide=0, short=0,
                spectral=0, fold_odd=0)
metric_count = 0

for N in range(3, max(args.center_max, args.metric_max)+1):
    for e in range(1, N+1, 2):
        a, d = N-e, N-2*e
        c = row(a, e)
        at = lambda k: c[k] if 0 <= k <= N else 0
        D = [at(k)**2-at(k-1)*at(k+1) for k in range(N+1)]
        E = [0]
        for v in D:
            E.append(E[-1]+v)

        if N <= args.center_max:
            for C in range(3, N+1):
                if (N-C) % 2:
                    continue
                x, n = (N-C)//2, C+1
                y = x+C
                S = E[y+1]-E[x]
                T = at(x)*at(y)-at(x-1)*at(y+1)
                assert S >= n*D[x] and D[x] == D[y] and S >= T
                if T <= 0:
                    branches["nonpositive"] += 1
                    continue
                assert at(x-1)**2 > at(x)**2
                if C*C >= N+4:
                    assert n*D[x] >= T
                    branches["wide"] += 1
                    continue
                f = c if d >= 0 else [
                    v if k % 2 == 0 else -v for k, v in enumerate(c)]
                if d >= 0 or C % 2 == 0:
                    if not long_sweep(f, x, y):
                        branches["short"] += 1
                        continue
                    b, key = domain(f, x, y+1), "spectral"
                else:
                    assert all(f[N-k] == f[k] for k in range(N+1))
                    b, key = domain(f, x-1, y+1), "fold_odd"
                assert b is not None
                L = b[1]-b[0]+1
                assert L <= n
                assert d*d*(L+1)**2 <= (N+1)**2*((L+1)**2-4)
                assert (n*n-1)*((N+3)**2-n*n) > n*n*d*d
                assert n*D[x] >= T
                branches[key] += 1

        if N <= args.metric_max:
            for C in range(3, N+1):
                for x in range((N-C+1)//2, N-C+1):
                    ok, L, eta, V, R, r, K = metric(N, d, x, C)
                    if not ok:
                        continue
                    y = x+C
                    S = E[y+1]-E[x]
                    T = at(x)*at(y)-at(x-1)*at(y+1)
                    assert S >= L*D[x]+(C+1-L)*D[y]
                    assert S*S >= eta*D[x]*D[y] >= T*T
                    b = T*T*(4*V+R*r)-K*D[x]*D[y]
                    assert b <= 0 or b*b <= 4*T**4*(R+r)**2*V
                    metric_count += 1

print("centered", sum(branches.values()), branches)
print("metric", metric_count)

def mul(p, q):
    out = {}
    for (i, j), b in p.items():
        for (k, l), d in q.items():
            out[i+k, j+l] = out.get((i+k, j+l), 0)+b*d
    return {k: v for k, v in out.items() if v}

def hp(k):
    return {(k-2*b-j, j): (-1)**b*comb(k+1-b, b)
            for b in range(k//2+1) for j in range(k+1-2*b)}

@lru_cache(None)
def cat(k):
    return 0 if k % 2 else comb(k, k//2)//(k//2+1)

@lru_cache(None)
def moment(r, i, j):
    return sum((-1)**k*comb(2*r, k)*cat(i+2*r-k)*cat(j+k)
               for k in range(2*r+1))

def direct(r, a, parts):
    p = {(0, 0): 1}
    for k in parts:
        p = mul(p, hp(k))
    p = mul(p, {(a-j, j): comb(a, j) for j in range(a+1)})
    value = sum(b*moment(r, i, j) for (i, j), b in p.items())
    assert value % 2 == 0
    return value//2

for r, a, parts in [
    (3, 8, (8, 8, 4)), (5, 16, (14, 14, 4)),
    (5, 2, (7, 7, 2)), (13, 3, (15, 14, 2))
]:
    u, v, w = parts
    e, C = 2*r-3, w+1
    N = a+e
    x = (N+u-v-C)//2
    y = x+C
    assert x+v+1 >= N+1
    c = row(a, e)
    at = lambda k: c[k] if 0 <= k <= N else 0
    S = sum(at(k)**2-at(k-1)*at(k+1) for k in range(x, y+1))
    T = at(x)*at(y)-at(x-1)*at(y+1)
    value = direct(r, a, parts)
    assert value == S-T
    print("Catalan", r, a, parts, value, S, T)
