import os, re, subprocess, time
from collections import defaultdict
from functools import lru_cache

# Census binary: build it from fm39/sec166_regime_split_census.cpp
#   g++ -O2 -std=c++17 -fopenmp -o flipx4 sec166_regime_split_census.cpp
# and pass its path in FLIPX4_BIN.
BIN = os.environ.get("FLIPX4_BIN", "./flipx4")

@lru_cache(None)
def profile(B):
    d = {(0, 0): 1}
    for z in B:
        n, eps = abs(z), (1 if z > 0 else -1)
        nd = defaultdict(int)
        for (a, b), v in d.items():
            for c in range(abs(a - n), a + n + 1, 2):
                nd[c, b] += v
            for c in range(abs(b - n), b + n + 1, 2):
                nd[a, c] += eps * v
        d = {key: value for key, value in nd.items() if value}
    return d

def canon(B):
    return tuple(sorted(B, key=lambda z: (abs(z), z)))

def fuses(a, k, p):
    return abs(a - k) <= p <= a + k and (a + k - p) % 2 == 0

def gp(B, p):
    return profile(canon(B)).get((p, 0), 0)

def child_value(C, k, p):
    # g_p(C S_k) = c_(p,k) + sum_(a: p in CG(a,k)) c_(a,0)
    d = profile(canon(C))
    return d.get((p, k), 0) + sum(
        v for (a, b), v in d.items() if b == 0 and fuses(a, k, p)
    )

def same_sign_pairs(B):
    labels = sorted(set(B), key=lambda z: (abs(z), z))
    out = []
    for i, x in enumerate(labels):
        if B.count(x) >= 2:
            out.append((x, x))
        for y in labels[i + 1:]:
            if (x > 0) == (y > 0):
                out.append((x, y))
    return out

def best_fusion(B, p):
    parent = gp(B, p)
    best = None
    for x, y in same_sign_pairs(B):
        C = list(B)
        C.remove(x)
        C.remove(y)
        a, b = abs(x), abs(y)
        for k in range(abs(a - b), a + b + 1, 2):
            if k == 0:
                continue
            candidate = (child_value(C, k, p), x, y, k)
            if best is None or candidate < best:
                best = candidate
    return parent, best

pat = re.compile(r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?) phi=(-?\d+)")
env = os.environ.copy()
env["OMP_NUM_THREADS"] = "4"
proc = subprocess.Popen([BIN, "36", "5"], stdout=subprocess.PIPE,
                        stderr=None, text=True, bufsize=1, env=env)

seen = passes = failures = 0
min_drop = None
start = time.monotonic()
for line in proc.stdout:
    if not line.startswith("NOFLIP "):
        continue
    m = pat.match(line)
    assert m, line
    W = int(m.group(1))
    p = abs(int(m.group(2)))
    B = canon(tuple(map(int, m.group(3).split())))
    phi = int(m.group(4))
    parent, best = best_fusion(B, p)
    seen += 1
    assert 2 * parent == phi
    if best is None or best[0] > parent:
        failures += 1
    else:
        passes += 1
        drop = parent - best[0]
        if min_drop is None or drop * min_drop[1] < min_drop[0] * parent:
            min_drop = (drop, parent, W, p, B, best)
    if seen % 200 == 0:
        print("progress", time.strftime("%H:%M:%S"),
              seen, "pass", passes, "fail", failures,
              "elapsed", round(time.monotonic() - start, 1), flush=True)

assert proc.wait() == 0
assert seen == 1904 and failures == 0
assert min_drop[:2] == (352, 429)
print("RESULT", {"records": seen, "pass": passes, "fail": failures,
                "minimum relative drop": min_drop,
                "seconds": round(time.monotonic() - start, 1)})
