cd /home/yang/q3adjoint
python3 -u - <<'FM187LAURENT'
from collections import defaultdict
from fractions import Fraction
from datetime import datetime, timezone

def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def wt(w):
    return sum(abs(z) for z in w)

def canon(w):
    return tuple(sorted(w, key=lambda z: (-abs(z), z)))

def interior_split(C):
    A, B = [], []
    wa = wb = 0
    for z in canon(C):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    return canon(A), canon(B)

def remove_pair(w, u, v):
    C = []
    ru = rv = False
    for z in w:
        if not ru and z == u:
            ru = True
        elif not rv and z == v:
            rv = True
        else:
            C.append(z)
    assert ru and rv
    return C

def laurent_table(w):
    p = {(0, 0): 1}
    for z in w:
        n = abs(z)
        eps = 1 if z > 0 else -1
        opts = [(i, 0, 1) for i in range(-n, n + 1, 2)]
        opts += [(0, j, eps) for j in range(-n, n + 1, 2)]
        q = defaultdict(int)
        for (r, s), c in p.items():
            for i, j, e in opts:
                q[(r + i, s + j)] += c * e
        p = {k: v for k, v in q.items() if v}
    cap = wt(w)
    f = {}
    for r in range(cap + 1):
        for s in range(cap + 1 - r):
            x = (p.get((r, s), 0) - p.get((r + 2, s), 0)
                 - p.get((r, s + 2), 0) + p.get((r + 2, s + 2), 0))
            if x:
                f[(r, s)] = x
    return f

tab_cache = {}
def tab(w):
    w = canon(w)
    if w not in tab_cache:
        tab_cache[w] = laurent_table(w)
    return tab_cache[w]

def pair_profile(w, u, v):
    C = remove_pair(w, u, v)
    A, B = interior_split(C)
    Bp = canon(list(B) + [u, v])
    FA, FB = tab(A), tab(Bp)
    cap = min(wt(A), wt(Bp))
    by_height = defaultdict(int)
    for (r, s), x in FA.items():
        if r + s <= cap:
            by_height[r + s] += x * FB.get((r, s), 0)
    parity = wt(A) % 2
    pref = 0
    mass = 0
    best = None
    best_data = None
    for T in range(parity, cap + 1, 2):
        layer = by_height[T]
        pref += layer
        mass += max(layer, 0)
        if mass:
            q = Fraction(pref, mass)
            if best is None or q < best:
                best, best_data = q, (T, pref, mass)
        elif pref < 0:
            best, best_data = None, (T, pref, mass)
            break
    return best, best_data

def all_pairs_best(w):
    types = list(dict.fromkeys(canon(w)))
    vals = []
    for i, u in enumerate(types):
        for v in types[i:]:
            if u == v and w.count(u) < 2:
                continue
            vals.append((pair_profile(w, u, v), (u, v)))
    return max(vals, key=lambda z: (Fraction(-10**100)
                                   if z[0][0] is None else z[0][0]))

seeds = [
    ("LP1", [-3]*14 + [2]*3 + [1]*2),
    ("LP2", [-3]*10 + [2]*7 + [1]*2),
]
for name, w in seeds:
    score, pair = all_pairs_best(w)
    print(stamp(), name, "best_pair", pair, "raw_min", score[0],
          "at(T,prefix,mass)", score[1])

w = [-3]*28 + [2]*9 + [1]*2
score, data = pair_profile(w, -3, -3)
print(stamp(), "GLOBAL_RAW_WORD", "W", wt(w), "n", len(w),
      "pair=(-3,-3)", "raw_min", score,
      "at(T,prefix,mass)", data)
assert score == Fraction(352474997782149202340448,
                        411061326807829926522413)
assert data == (49, 704949995564298404680896,
                822122653615659853044826)

def reflection_variants(w):
    r = [-z for z in w]
    out = []
    for label in sorted(set(abs(z) for z in r)):
        idx = next((i for i,z in enumerate(r)
                    if abs(z)==label and z<0), None)
        if idx is None:
            idx = next(i for i,z in enumerate(r)
                       if abs(z)==label and z>0)
        q = r.copy()
        q[idx] = -q[idx]
        out.append((label, q))
    return out

for name, w in seeds:
    for label, rw in reflection_variants(w):
        types = list(dict.fromkeys(canon(rw)))
        records = []
        for i,u in enumerate(types):
            for v in types[i:]:
                if u == v and rw.count(u) < 2:
                    continue
                q, data = pair_profile(rw, u, v)
                fail = (q is not None and q < 0) or (
                    q is None and data is not None and data[1] < 0)
                records.append((q, data, (u,v), fail))
        good = sum(not rec[3] for rec in records)
        best = max(records, key=lambda rec:
                   Fraction(-10**100) if rec[0] is None else rec[0])
        print(stamp(), name+"_reflection_flip_label_"+str(label),
              "minus_count",sum(z<0 for z in rw),
              "best_pair",best[2],"best_raw",best[0],
              "at",best[1],"passing_pairs",good,"of",len(records))
FM187LAURENT
