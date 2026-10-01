import random, sys, time
from collections import defaultdict
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
def score(bg, ins):
    # bg: all-minus plus <= 1 odd plus; ins: two minus labels; return min relative increase over admissible p
    pos = [v for v in bg if v > 0]
    if len(pos) > 1 or any(v % 2 == 0 for v in pos): return None
    if any((-v in bg) for v in pos) or any(v in pos for v in [-i for i in ins]): return None
    g0 = gv(bg); g1 = step(step(g0, abs(ins[0]), -1), abs(ins[1]), -1)
    mx = max([abs(v) for v in bg + ins] + [3]); W1 = sum(abs(v) for v in bg + ins)
    best = None
    for p in range(mx, W1+1):
        if (W1 - p) % 2: continue
        a = g0.get((p, 0), 0); b = g1.get((p, 0), 0)
        if b < a: return (-1.0, p, a, b)
        if a > 0:
            r = (b - a)/a
            if best is None or r < best[0]: best = (r, p, a, b)
    return best
def run(seed):
    rng = random.Random(seed); t0 = time.time(); bestg = None
    while time.time() - t0 < float(sys.argv[1]):
        bg = [-rng.randint(1, 12) for _ in range(rng.randint(1, 8))]
        if rng.random() < 0.5:
            o = rng.choice([1, 3, 5, 7, 9]); 
            if -o not in bg: bg.append(o)
        ins = [-rng.randint(1, 12), -rng.randint(1, 12)]
        cur = score(bg, ins)
        if cur is None: continue
        for it in range(120):
            nb, ni = list(bg), list(ins); op = rng.random()
            if op < 0.3 and len(nb) > 1: nb.pop(rng.randrange(len(nb)))
            elif op < 0.6: nb.append(-rng.randint(1, 14))
            elif op < 0.8: ni[rng.randrange(2)] = -rng.randint(1, 14)
            else:
                i = rng.randrange(len(nb)); nb[i] = (1 if nb[i] > 0 else -1)*max(1, abs(nb[i]) + rng.choice([-1, 1]))
            if sum(abs(v) for v in nb + ni) > 110: continue
            sc = score(nb, ni)
            if sc is None: continue
            if sc[0] <= cur[0]: bg, ins, cur = nb, ni, sc
            if cur[0] < 0: return (cur, sorted(bg), ins)
        if bestg is None or cur[0] < bestg[0][0]: bestg = (cur, sorted(bg), ins)
    return bestg
if __name__ == "__main__":
    with Pool(32) as pool: res = [r for r in pool.map(run, range(100, 132)) if r]
    res.sort(key=lambda r: r[0][0])
    for r in res[:6]: print(f"min rel increase {r[0][0]:.4f} at p {r[0][1]}  bg {r[1]}  insert {r[2]}")
    print("VIOLATION FOUND" if res and res[0][0][0] < 0 else "no violation")
