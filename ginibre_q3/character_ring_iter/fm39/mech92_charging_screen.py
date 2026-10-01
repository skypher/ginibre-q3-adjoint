# Plausibility test: walk expansion F_B(k) = sum M_B(k;i,j) W(i,j) and a crude (E)-charging bound.
from math import comb
import random, sys

def row(a, e):
    c = [1]
    for _ in range(a): c = [x + y for x, y in zip(c + [0], [0] + c)]
    for _ in range(e): c = [x - y for x, y in zip(c + [0], [0] + c)]
    return c

def get(c, j): return c[j] if 0 <= j < len(c) else 0

def Kshift(c):  # (Kf)_j = f_(j-1) + f_(j+1), on indices -B..N+B -> use dict
    lo, hi = min(c), max(c)
    return {j: c.get(j-1, 0) + c.get(j+1, 0) for j in range(lo-1, hi+2)}

def F_direct(a, e, B, k):
    c0 = {j: v for j, v in enumerate(row(a, e))}
    vs = [c0]
    for h in range(B+1): vs.append(Kshift(vs[-1]))
    def D(h, k): v = vs[h]; return v.get(k, 0)**2 - v.get(k-1, 0)*v.get(k+1, 0)
    def T(b, k): return sum(comb(b, h) * 2**(b-h) * D(h, k) for h in range(b+1))
    return T(B, k) - T(B, k+1)

def walks(k, B):
    M = {(k, k+1): 1}
    for _ in range(B):
        nxt = {}
        for (i, j), w in M.items():
            g = j - i
            moves = [((i-1, j-1), 1), ((i-1, j+1), 1), ((i+1, j+1), 1), ((i, j), 1 if g == 1 else 2)]
            if g >= 3: moves.append(((i+1, j-1), 1))
            for key, wt in moves: nxt[key] = nxt.get(key, 0) + w*wt
        M = nxt
    return M

def W(c, i, j):
    pi = (get(c, i), get(c, i-1) + get(c, i+1)); pj = (get(c, j), get(c, j-1) + get(c, j+1))
    return pi[0]*pj[1] - pi[1]*pj[0]

def Dt(c, l): return get(c, l)**2 - get(c, l-1)*get(c, l+1)

random.seed(1)
bad_expansion = 0; tests = 0; charge_ok = 0; charge_fail = 0; negs = 0
for trial in range(4000):
    a = random.randint(2, 40); e = random.randint(2, 40); N = a + e
    B = random.randint(1, 8)
    # consumer k = (N+p)/2 with p >= 1, N+p even; restrict to right-half windows
    p = random.randrange(1 if N % 2 else 2, N + 2*B + 2, 2)
    k = (N + p)//2
    c = row(a, e)
    M = walks(k, B)
    Fw = sum(w * W(c, i, j) for (i, j), w in M.items() if i < j)
    Fd = F_direct(a, e, B, k)
    tests += 1
    if Fw != Fd: bad_expansion += 1; continue
    if Fd < 0: negs += 1
    # crude charge: adjacent chords exact; longer chords: max(W,0) dropped, negative bounded by -(D_i - D_j)
    lowest = min(i for (i, j) in M); 
    if 2*lowest < N: continue  # charging uses (E) only on the right half
    L = 0
    for (i, j), w in M.items():
        if j - i == 1: L += w * (Dt(c, i) - Dt(c, j))
        else: L -= w * max(0, Dt(c, i) - Dt(c, j))
    if L >= 0: charge_ok += 1
    else: charge_fail += 1
print("tests", tests, "expansion mismatches", bad_expansion, "negative F", negs)
print("right-half windows: crude charge >= 0:", charge_ok, " < 0:", charge_fail)
