# Adversarial test of Venus's duplicate-first rule (FM-SEC148) on pair-free backgrounds, W up to WMAX.
import random, sys, time
from collections import defaultdict, Counter
from multiprocessing import Pool
def step(st, n, e):
    nx = defaultdict(int)
    for (s, t), v in st.items():
        for s2 in range(abs(s-n), s+n+1, 2): nx[s2, t] += v
        for t2 in range(abs(t-n), t+n+1, 2): nx[s, t2] += e*v
    return nx
def gv(lab):
    st = {(0, 0): 1}
    for x in lab: st = step(st, abs(x), 1 if x > 0 else -1)
    return st
def rule1(B):
    c = Counter(B)                                   # signed labels; pair-free so value -> one sign
    evp = sorted(v for v in c if v > 0 and v % 2 == 0)
    if MODE in ("evenplus_first", "partb") and evp: return [evp[0]]
    best = None
    for v, f in c.items():
        if f >= 2 and (best is None or f > best[0] or (f == best[0] and abs(v) < abs(best[1]))): best = (f, v)
    if best: return [best[1], best[1]]
    vals = sorted(c, key=abs)
    pairs = [(a, b) for i, a in enumerate(vals) for b in vals[i+1:] if (abs(a) + abs(b)) % 2 == 0]
    if pairs:
        a, b = max(pairs, key=lambda ab: (max(abs(ab[0]), abs(ab[1])), min(abs(ab[0]), abs(ab[1]))))
        return [a, b]
    ev = [v for v in vals if abs(v) % 2 == 0]
    return [ev[0]] if ev else None
def margin(B):
    R = rule1(B)
    if R is None: return None
    sm = list(B)
    for r in R: sm.remove(r)
    g1 = gv(B); g0 = gv(sm)
    mx = max([abs(v) for v in B] + [3]); W = sum(abs(v) for v in B); best = None
    for p in range(mx, W + 1):
        if (W - p) % 2: continue
        a = g0.get((p, 0), 0); b = g1.get((p, 0), 0)
        if b < a: return (-1.0, p, R, a, b)
        r = (b - a) / max(1, a)
        if best is None or r < best[0]: best = (r, p, R, a, b)
    return best
def pf(B): return not any(-v in B for v in B)
WMAX = int(sys.argv[2])
MODE = sys.argv[3] if len(sys.argv) > 3 else 'plain'
def run(seed):
    rng = random.Random(seed); t0 = time.time(); bestg = None
    while time.time() - t0 < float(sys.argv[1]):
        sg = {}; B = []
        mult = rng.randint(1, 4)
        for v in rng.sample(range(1, 16), rng.randint(1, 7)):
            sg[v] = rng.choice([1, -1])
            if sg[v] > 0 and v % 2 == 0: sg[v] = -1
            B += [sg[v]*v]*max(1, mult + rng.choice([-1, 0, 0, 1]))
        if len(B) < 2 or sum(abs(v) for v in B) > WMAX: continue
        cur = margin(B)
        if cur is None: continue
        for it in range(150):
            nb = list(B); op = rng.random()
            if op < 0.2 and len(nb) > 2: nb.pop(rng.randrange(len(nb)))
            elif op < 0.45: nb.append(rng.choice(nb))
            elif op < 0.6:
                v = rng.randint(1, 18); s = rng.choice([1, -1])
                if -s*v in nb: s = -s
                nb.append(s*v)
            elif op < 0.75:
                v = abs(rng.choice(nb)); nb = [(-x if abs(x) == v else x) for x in nb]
            else:
                i = rng.randrange(len(nb)); x = nb[i]; y = (1 if x > 0 else -1)*max(1, abs(x) + rng.choice([-1, 1]))
                if -y in nb: continue
                nb[i] = y
            if not pf(nb) or len(nb) < 2 or sum(abs(v) for v in nb) > WMAX: continue
            if MODE == "partb" and any(v > 0 and v % 2 == 0 for v in nb): continue
            sc = margin(nb)
            if sc is None: continue
            if sc[0] <= cur[0]: B, cur = nb, sc
            if cur[0] < 0: return (cur, sorted(B))
        if bestg is None or cur[0] < bestg[0][0]: bestg = (cur, sorted(B))
    return bestg
if __name__ == "__main__":
    with Pool(32) as pool: res = [r for r in pool.map(run, range(800, 832)) if r]
    res.sort(key=lambda r: r[0][0])
    for r in res[:6]: print(f"min rel increase {r[0][0]:.4f} at p {r[0][1]} removal {r[0][2]} background {r[1]}")
    print(f"[{MODE}]", "VIOLATION FOUND" if res and res[0][0][0] < 0 else "no violation")
