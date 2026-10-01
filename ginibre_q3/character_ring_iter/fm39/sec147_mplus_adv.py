# Adversarial falsification of (M+): pair-free background B, insertion = one even plus label or two odd plus labels.
import random, sys, time
from collections import defaultdict
from multiprocessing import Pool
exec(open(__file__.replace('sec147_mplus_adv.py', 'sec147_mminus_adv.py')).read().split('def score')[0])
MODE = sys.argv[2] if len(sys.argv) > 2 else "all"
def pairfree(l): return not any((-v in l) for v in l)
def score(bg, ins):
    big = bg + ins
    if not pairfree(big): return None
    if MODE != "E1" and any(v < 0 for v in ins): return None
    if len(ins) == 1 and ins[0] % 2: return None
    if len(ins) == 2 and (ins[0] % 2 == 0 or ins[1] % 2 == 0): return None
    if MODE in ("even", "E1") and len(ins) != 1: return None
    if MODE == "oddlargest":
        if len(ins) != 2: return None
        oddplus = [v for v in bg if v > 0 and v % 2 == 1]
        if oddplus and min(ins) < max(oddplus): return None
    g0 = gv(bg); g1 = gv(big)
    mx = max([abs(v) for v in big] + [3]); W = sum(abs(v) for v in big)
    best = None
    for p in range(mx, W+1):
        if (W - p) % 2: continue
        a = g0.get((p, 0), 0); b = g1.get((p, 0), 0)
        if b < a: return (-1.0, p, a, b)
        if a > 0:
            r = (b - a)/a
            if best is None or r < best[0]: best = (r, p, a, b)
    return best
def rand_ins(rng, bg):
    for _ in range(30):
        if MODE == "E1": ins = [rng.choice([1, -1])*2*rng.randint(1, 7)]
        elif MODE == "even" or (MODE != "oddlargest" and rng.random() < 0.5): ins = [2*rng.randint(1, 7)]
        elif MODE == "oddlargest":
            mo = max([v for v in bg if v > 0 and v % 2 == 1] + [1])
            ins = [mo + 2*rng.randint(0, 4), mo + 2*rng.randint(0, 4)]
        else: ins = [2*rng.randint(0, 6)+1, 2*rng.randint(0, 6)+1]
        if pairfree(bg + ins): return ins
    return None
def run(seed):
    rng = random.Random(seed); t0 = time.time(); bestg = None
    while time.time() - t0 < float(sys.argv[1]):
        sg = {}; bg = []
        for v in rng.sample(range(1, 13), rng.randint(1, 6)):
            sg[v] = rng.choice([1, -1]); bg += [sg[v]*v]*rng.randint(1, 4)
        ins = rand_ins(rng, bg)
        if ins is None: continue
        cur = score(bg, ins)
        if cur is None: continue
        for it in range(120):
            nb, ni = list(bg), list(ins); op = rng.random()
            if op < 0.25 and len(nb) > 1: nb.pop(rng.randrange(len(nb)))
            elif op < 0.5: nb.append(rng.choice(nb))
            elif op < 0.65:
                v = rng.randint(1, 14); s = rng.choice([1, -1]); nb.append(s*v)
            elif op < 0.8:
                v = abs(rng.choice(nb)); nb = [(-x if abs(x) == v else x) for x in nb]
            else:
                ni = rand_ins(rng, nb) or ni
            if sum(abs(v) for v in nb + ni) > 110: continue
            sc = score(nb, ni)
            if sc is None: continue
            if sc[0] <= cur[0]: bg, ins, cur = nb, ni, sc
            if cur[0] < 0: return (cur, sorted(bg), ins)
        if bestg is None or cur[0] < bestg[0][0]: bestg = (cur, sorted(bg), ins)
    return bestg
if __name__ == "__main__":
    with Pool(32) as pool: res = [r for r in pool.map(run, range(300, 332)) if r]
    res.sort(key=lambda r: r[0][0])
    for r in res[:6]: print(f"min rel increase {r[0][0]:.4f} at p {r[0][1]}  bg {r[1]}  insert {r[2]}")
    print(f"[{MODE}]", "VIOLATION FOUND" if res and res[0][0][0] < 0 else "no violation")
