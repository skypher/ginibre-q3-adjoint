from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import time

def stamp():
    return time.strftime("%H:%M:%S", time.gmtime())

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def mul_raw(A, B):
    out = defaultdict(int)
    for (r, s), va in A.items():
        for (u, v), vb in B.items():
            for x in cg(r, u):
                for y in cg(s, v):
                    out[x, y] += va * vb
    return {k: value for k, value in out.items() if value}

def multiply_d(tab):
    out = defaultdict(int)
    for (r, s), value in tab.items():
        for x in cg(r, 1):
            out[x, s] += value
        for y in cg(s, 1):
            out[r, y] -= value
    return {k: value for k, value in out.items() if value}

def mul_factor(tab, z):
    n = abs(z)
    eps = 1 if z > 0 else -1
    out = defaultdict(int)
    for (r, s), value in tab.items():
        for t in cg(r, n):
            out[t, s] += value
        for t in cg(s, n):
            out[r, t] += eps * value
    return {k: v for k, v in out.items() if v}

@lru_cache(None)
def table_cached(word):
    tab = {(0, 0): 1}
    for z in word:
        tab = mul_factor(tab, z)
    return tab

def table(word):
    return table_cached(tuple(sorted(word, key=lambda z: (-abs(z), z))))

def quotient_table(word):
    q = {(0, 0): 1}
    for z in word:
        n = abs(z)
        factor = ({(j, n-1-j): 1 for j in range(n)} if z < 0
                  else {(n, 0): 1, (0, n): 1})
        q = mul_raw(q, factor)
    return q

@lru_cache(None)
def character_coefficients(word, extra_d):
    q = quotient_table(word)
    for _ in range(extra_d):
        q = multiply_d(q)
    alternating = multiply_d(q)
    return tuple(sorted(((i-1, j), value)
                        for (i, j), value in alternating.items() if i > j))

def interior_cut(word):
    A, B = [], []
    wa = wb = 0
    for z in sorted(word, key=lambda z: (-abs(z), z)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    return A, B

def pair_prefix(L, i, j):
    C = [z for k, z in enumerate(L) if k not in (i, j)]
    A, B = interior_cut(C)
    Bp = B + [L[i], L[j]]
    FA, FB = table(A), table(Bp)
    layer = defaultdict(int)
    for (r, s), value in FA.items():
        layer[r+s] += value * FB.get((r, s), 0)
    running = 0
    minimum = None
    argmin = None
    for height in sorted(layer):
        running += layer[height]
        if minimum is None or running < minimum:
            minimum, argmin = running, height
        assert running >= 0, (L, i, j, height, running)
    imbalance = abs(sum(map(abs, A)) - sum(map(abs, Bp)))
    return minimum, argmin, running, A, Bp, imbalance

def audit(name, L):
    start = time.time()
    count = 0
    least = None
    where = None
    for i, j in combinations(range(len(L)), 2):
        minimum, height, phi, _, _, _ = pair_prefix(L, i, j)
        count += 1
        if least is None or minimum < least:
            least = minimum
            where = (L[i], L[j], height, phi)
    print(stamp(), name, "pairs", count, "least_prefix", least,
          "at", where, "seconds", round(time.time()-start, 2), flush=True)
    return count

def top_pair(B):
    same = []
    all_pairs = []
    for i, j in combinations(range(len(B)), 2):
        key = (abs(B[i]) + abs(B[j]), max(abs(B[i]), abs(B[j])))
        all_pairs.append((key, i, j))
        if (abs(B[i]) - abs(B[j])) % 2 == 0:
            same.append((key, i, j))
    options = same if same else all_pairs
    assert options
    _, i, j = max(options)
    return i, j

def l_image(lam, T):
    alpha, beta = lam
    ell = alpha - beta
    if T < ell:
        return None
    j = min(beta, (T-ell)//2)
    return (ell+j, j)

def gamma_value(ca, cb, T, mode):
    if mode == "11":
        return sum(a * cb.get(lam, 0) for lam, a in ca.items()
                   if sum(lam) <= T-1)
    if mode == "02":
        return sum(a * cb.get(l_image(lam, T), 0)
                   for lam, a in ca.items()
                   if l_image(lam, T) is not None)
    if mode == "20":
        return sum(b * ca.get(l_image(lam, T), 0)
                   for lam, b in cb.items()
                   if l_image(lam, T) is not None)
    raise ValueError(mode)

def check_channel_domination(name, L):
    i, j = top_pair(L[:-1])
    A, B = interior_cut([z for k, z in enumerate(L) if k not in (i, j)])
    Bp = B + [L[i], L[j]]
    a, b = sum(z < 0 for z in A), sum(z < 0 for z in Bp)
    options = []
    if a % 2 == 1:
        options.append((1, 1, "11"))
    elif a == b == 0:
        print(stamp(), name, "all-plus coefficient case", flush=True)
        return
    else:
        if b >= 2:
            options.append((0, 2, "02"))
        if a >= 2:
            options.append((2, 0, "20"))
    assert options and a % 2 == b % 2
    results = []
    total = sum(map(abs, L))
    fa, fb = table(A), table(Bp)
    for ra, rb, mode in options:
        ca = dict(character_coefficients(tuple(A), a-ra))
        cb = dict(character_coefficients(tuple(Bp), b-rb))
        minimum = None
        for T in range(total+3):
            gamma = gamma_value(ca, cb, T, mode)
            prefix = sum(x*fb.get(k, 0) for k, x in fa.items() if sum(k) <= T)
            assert 2*gamma == prefix, (name, mode, T, 2*gamma, prefix)
            if minimum is None or gamma < minimum:
                minimum = gamma
        assert minimum >= 0, (name, mode, minimum)
        results.append((mode, (ra, rb), minimum))
    print(stamp(), name, "pair", (L[i], L[j]), "minus", (a, b),
          "channel-domination", results, flush=True)

def grouped_coefficients(coeffs):
    grouped = defaultdict(dict)
    for (alpha, beta), value in coeffs.items():
        q = alpha + beta
        grouped[q][beta] = value
    return {q: tuple(grouped[q].get(beta, 0)
                     for beta in range(q//2+1))
            for q in sorted(grouped)}

# Exact segregated split and its L_T channel domination failure.
bad_A = [1, -2, -4, -4, -4]
bad_B = [1] + [3]*6 + [5]*4 + [6]
TA, TB = table(bad_A), table(bad_B)
bad_layers = defaultdict(int)
for (r, s), value in TA.items():
    bad_layers[r+s] += value * TB.get((r, s), 0)
pi11 = sum(value for h, value in bad_layers.items() if h <= 11)
phi_bad = sum(bad_layers.values())
assert (pi11, phi_bad) == (-24695910, 120541550)
ca_bad = dict(character_coefficients(tuple(bad_A), 2))
cb_bad = dict(character_coefficients(tuple(bad_B), 0))
positive = negative = 0
for lam, bcoef in cb_bad.items():
    image = l_image(lam, 11)
    if image is not None:
        product = bcoef * ca_bad.get(image, 0)
        if product >= 0:
            positive += product
        else:
            negative -= product
assert (positive, negative, 2*(positive-negative)) == (
    288054364, 300402319, -24695910)
print(stamp(), "segregated T=11 channel masses", positive, negative,
      "Pi_11", 2*(positive-negative), flush=True)

# Full W1 character decompositions and the J_(T-1) channel payment.
A1 = [6, -4, -4, 3, 3, 3, -2]
B1 = [5, 5, -4, 3, 3, 3, 1, 1, 5, 5]
GA = dict(character_coefficients(tuple(A1), 2))
GB = dict(character_coefficients(tuple(B1), 0))
vA, vB = grouped_coefficients(GA), grouped_coefficients(GB)
assert (len(GA), sum(x < 0 for x in GA.values()),
        len(GB), sum(x < 0 for x in GB.values())) == (89, 53, 170, 1)
assert max(a+b for a,b in GA) == 24 and max(a+b for a,b in GB) == 34
assert GB[(16,16)] == -12 and (16,16) not in GA
print(stamp(), "W1 G_A vectors beta=0..floor(q/2)", vA, flush=True)
print(stamp(), "W1 G_B vectors beta=0..floor(q/2)", vB, flush=True)

stats = []
cumulative = 0
for q in sorted(set(a+b for a,b in GA) | set(a+b for a,b in GB)):
    products = [GA[lam]*GB[lam] for lam in GA.keys() & GB.keys()
                if sum(lam) == q]
    pos = sum(x for x in products if x > 0)
    neg = -sum(x for x in products if x < 0)
    net = pos-neg
    cumulative += net
    if products:
        stats.append((q, pos, neg, net, cumulative))
expected = [
 (0,1969016,0,1969016,1969016),
 (2,7037150,796404,6240746,8209762),
 (4,17421695,3649268,13772427,21982189),
 (6,28023034,12246718,15776316,37758505),
 (8,26987352,13008961,13978391,51736896),
 (10,37147721,29489757,7657964,59394860),
 (12,32745706,28963604,3782102,63176962),
 (14,15356598,15817962,-461364,62715598),
 (16,13882214,15084861,-1202647,61512951),
 (18,6269568,6980633,-711065,60801886),
 (20,1331711,1741321,-409610,60392276),
 (22,714640,829336,-114696,60277580),
 (24,147848,154653,-6805,60270775),
]
assert stats == expected, stats
assert 2*min(x[4] for x in stats) == 3938032
assert 2*stats[-1][4] == 120541550
print(stamp(), "W1 matched-channel (q,positive,negative,net,cumulative)", stats,
      flush=True)

# Three SEC181 witnesses, F1, the census tight case, and the runs.
words = [
 ("SEC181-W1", [1,1,-2]+[3]*6+[-4]*3+[5]*4+[6]),
 ("SEC181-W2", [1,1,-2]+[3]*5+[-4]*3+[5]*5+[6]),
 ("SEC181-W3", [1,1,1,-2,2]+[3]*5+[-4]*3+[4]*2+[5]*2+[8]),
]
checks = 0
for name, word in words:
    checks += audit(name, word)
F1 = [1,1,-2]+[3]*4+[-4]*3+[5]*4+[8]
checks += audit("F1", F1)
tight = [-1,-3,-3,-4,-5,-5,-6,-7]
checks += audit("census-tight", tight)

for name, word in words:
    check_channel_domination(name, word)
check_channel_domination("F1", F1)
check_channel_domination("census-tight", tight)

run_cases = 0
for k in range(2, 11):
    sign = (-1)**k
    for p0 in range(k+1, k+4):
        if (k*(k+1)//2+p0) % 2:
            continue
        run = [-n for n in range(1,k+1)] + [sign*p0]
        checks += audit(f"run-k{k}-p{sign*p0}", run)
        check_channel_domination(f"run-k{k}-p{sign*p0}", run)
        run_cases += 1

assert checks == 890 and run_cases == 14
print(stamp(), "exact pair-prefix checks", checks, "run cases", run_cases,
      "cache entries", table_cached.cache_info().currsize, "PASS", flush=True)
