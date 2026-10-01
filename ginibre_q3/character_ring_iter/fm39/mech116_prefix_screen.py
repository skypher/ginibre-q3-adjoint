import argparse
from functools import lru_cache
from itertools import combinations_with_replacement, product
from fractions import Fraction
from collections import defaultdict

parser = argparse.ArgumentParser(description="Exact FM3 toggle-prefix screens.")
parser.parse_args()

@lru_cache(None)
def fusion(ns):
    state = {0: 1}
    for n in ns:
        nxt = {}
        for a, v in state.items():
            for b in range(abs(a-n), a+n+1, 2):
                nxt[b] = nxt.get(b, 0) + v
        state = nxt
    return state

def inv(ns):
    return fusion(tuple(sorted(abs(n) for n in ns))).get(0, 0)

def split_products(ns):
    L = len(ns)
    out = []
    for s in range(1 << L):
        A = tuple(abs(ns[j]) for j in range(L) if s >> j & 1)
        B = tuple(abs(ns[j]) for j in range(L) if not (s >> j & 1))
        out.append(inv(A) * inv(B))
    return out

def toggle_layers(ns, eps, i, fs):
    L = len(ns)
    minus = sum(1 << j for j, e in enumerate(eps) if e < 0)
    T = [0] * L
    for s in range(1 << L):
        if s >> i & 1:
            continue
        parity = (s & minus).bit_count() & 1
        T[s.bit_count()] += (-1 if parity else 1) * (
            fs[s] - fs[s | (1 << i)]
        )
    pref = []
    z = 0
    for v in T:
        z += v
        pref.append(z)
    return T, pref

def full_split_sum(ns, eps, fs):
    minus = sum(1 << j for j, e in enumerate(eps) if e < 0)
    return sum((-1 if (s & minus).bit_count() & 1 else 1) * fs[s]
               for s in range(1 << len(ns)))

def sign_budget(ns, eps):
    C = Fraction(0)
    for j in range(len(ns)-1):
        c = sum(1 for s in range(1 << (j+1))
                if s >> j & 1 and
                sum(eps[k] < 0 for k in range(j+1) if s >> k & 1) % 2)
        C += Fraction(c, ns[j] + 1)
    return C

def hpoly(word, signs, a, b):
    # Coefficients of E[U_a(x) U_b(y) prod_j(U_nj(x)+eps_j z U_nj(y))].
    L = len(word)
    out = [0] * (L+1)
    for mask in range(1 << L):
        A = tuple(abs(word[j]) for j in range(L) if mask >> j & 1)
        B = tuple(abs(word[j]) for j in range(L) if not (mask >> j & 1))
        sgn = -1 if sum(signs[j] < 0 for j in range(L)
                        if mask >> j & 1) % 2 else 1
        out[mask.bit_count()] += sgn * fusion(B).get(a, 0) * fusion(A).get(b, 0)
    return out

def poly_add(a, b, scale=1, shift=0):
    out = [0] * max(len(a), len(b) + shift)
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i+shift] += scale*v
    return out

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def gvalue(word, p):
    state = {(0, 0): 1}
    for sn in word:
        n = abs(sn)
        eps = 1 if sn > 0 else -1
        nxt = defaultdict(int)
        for (a, b), v in state.items():
            for c in range(abs(a-n), a+n+1, 2):
                nxt[c, b] += v
            for d in range(abs(b-n), b+n+1, 2):
                nxt[a, d] += eps*v
        state = nxt
    return state.get((p, 0), 0)

def screen(values, maxL):
    profiles = with_minus = toggles = neg_layers = neg_prefixes = 0
    min_prefix = None
    for L in range(2, maxL+1):
        for ns in combinations_with_replacement(values, L):
            if not inv(ns):
                continue
            vals = tuple(sorted(set(ns)))
            fs = split_products(ns)
            for ss in product((1, -1), repeat=len(vals)):
                d = dict(zip(vals, ss))
                eps = tuple(d[n] for n in ns)
                if sum(e < 0 for e in eps) % 2:
                    continue
                profiles += 1
                if any(e < 0 for e in eps):
                    with_minus += 1
                for i, e in enumerate(eps):
                    if e > 0:
                        continue
                    toggles += 1
                    T, pref = toggle_layers(ns, eps, i, fs)
                    assert sum(T) == full_split_sum(ns, eps, fs)
                    neg_layers += any(v < 0 for v in T)
                    neg_prefixes += any(v < 0 for v in pref)
                    local = min(pref)
                    if min_prefix is None or local < min_prefix[0]:
                        min_prefix = (local, (ns, eps, i, T, pref))
    return profiles, with_minus, toggles, neg_layers, neg_prefixes, min_prefix

core = screen(range(3, 8), 7)
assert core[:5] == (1900, 1511, 5492, 1011, 0)
print("pure-core labels 3..7, L<=7:", core[:5], "min prefix:", core[5])

small = screen(range(1, 7), 8)
assert small[:5] == (10627, 9126, 37382, 13178, 0)
assert small[5][0] == 1
print("pair-free labels 1..6, L<=8:", small[:5], "min prefix:", small[5])

# Unsaturated pure-core delta=8 frontier with 2..6 background factors.
delta = 8
profiles = selected = failures = 0
min_target = None
for k in range(2, 7):
    for bg in combinations_with_replacement(range(3, 11), k):
        mx = max(bg)
        if mx > delta:
            continue
        p = sum(bg) - 2*delta
        if p < mx or p < 3:
            continue
        vals = tuple(sorted(set(bg)))
        for ss in product((1, -1), repeat=len(vals)):
            d = dict(zip(vals, ss))
            ebg = tuple(d[n] for n in bg)
            ep = 1 if sum(e < 0 for e in ebg) % 2 == 0 else -1
            if p in d and d[p] != ep:
                continue
            ns, eps = bg + (p,), ebg + (ep,)
            if not inv(ns):
                continue
            profiles += 1
            if ep < 0:
                i = k
            else:
                ids = [j for j, e in enumerate(ebg) if e < 0]
                if not ids:
                    continue
                i = ids[0]
            selected += 1
            fs = split_products(ns)
            T, pref = toggle_layers(ns, eps, i, fs)
            assert sum(T) == full_split_sum(ns, eps, fs)
            failures += any(v < 0 for v in pref)
            if min_target is None or min(pref) < min_target[0]:
                min_target = (min(pref), (bg, ebg, p, ep, i, T, pref))
assert (profiles, selected, failures) == (7678, 6948, 0)
assert min_target[0] == 9
print("delta=8 pure-core unsaturated screen:", profiles, selected,
      failures, "min prefix:", min_target)

# Exact C>1 frontier word: B=(-3,-3,-4,-6,-7), distinguished -7.
word = (3, 3, 4, 6, 7, 7)
signs = (-1,) * 6
i = 5
fs = split_products(word)
T, pref = toggle_layers(word, signs, i, fs)
assert inv(word) == 71 and sign_budget(word, signs) == Fraction(173, 70)
assert T == [71, 3, -1, -1, 3, 71]
assert pref == [71, 74, 73, 72, 75, 146]
assert full_split_sum(word, signs, fs) == 146
bg = (-3, -3, -4, -6, -7)
assert gvalue(bg, 7) == 73
assert fusion(tuple(abs(x) for x in bg)).get(7, 0) == 71
print("delta=8 frontier:", "C=", sign_budget(word, signs),
      "layers=", T, "prefixes=", pref, "g7=", gvalue(bg, 7),
      "m7=", fusion(tuple(abs(x) for x in bg)).get(7, 0))

# Exact generator action on H_R(a,b;z) for one repeated -7 insertion.
R = (-3, -3, -4, -6)
eR = (-1,) * 4
m, eta, n = 7, -1, 7
child = R + (m,)
echild = eR + (eta,)
q_direct = poly_add(hpoly(child, echild, n, 0),
                    hpoly(child, echild, 0, n), scale=-1)
q_update = [0] * (len(R)+2)
for c in cg(n, m):
    q_update = poly_add(q_update, hpoly(R, eR, c, 0))
q_update = poly_add(q_update, hpoly(R, eR, n, m), scale=eta, shift=1)
q_update = poly_add(q_update, hpoly(R, eR, m, n), scale=-1)
for c in cg(n, m):
    q_update = poly_add(q_update, hpoly(R, eR, 0, c), scale=-eta, shift=1)
assert q_direct == q_update == T
print("generator update exact:", q_update)

# Exact low-label extension obstruction to value monotonicity.
low_bg = (-9,) * 5 + (4,) * 6 + (12,) * 3
p = 13
before = gvalue(low_bg, p)
after = gvalue(low_bg + (1, 1), p)
assert (before, after, after-before) == (
    90245874567, 80642598148, -9603276419
)
print("adding +1,+1: g13 before, after, change:",
      before, after, after-before)
print("PASS")