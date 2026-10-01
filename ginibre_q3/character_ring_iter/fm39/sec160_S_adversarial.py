import random, subprocess, sys, time
from multiprocessing import Pool
EXE = '/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/dtest/ruleS1'
def score(B):
    out = subprocess.run([EXE] + [str(x) for x in B], capture_output=True, text=True).stdout.split()
    W, tot, neg, mn = int(out[0]), int(out[1]), int(out[2]), float(out[3])
    if tot == 0: return None
    return (-1.0 if neg else mn), neg
def pf(B): return not any(-v in B for v in B)
WMAX = int(sys.argv[2])
def run(seed):
    rng = random.Random(seed); t0 = time.time(); best = None
    while time.time() - t0 < float(sys.argv[1]):
        sg = {}; B = []; mu = rng.randint(1, 4); flat = rng.random() < 0.5
        for v in rng.sample(range(1, 26), rng.randint(3, 11)):
            sg[v] = rng.choice([1, -1]); B += [sg[v]*v]*(mu if flat else rng.randint(1, 6))
        if rng.random() < 0.5: B += [sg.setdefault(1, rng.choice([1, -1]))]*rng.randint(5, 40)
        if not pf(B) or sum(abs(x) for x in B) > WMAX: continue
        cur = score(B)
        if cur is None: continue
        for it in range(50):
            nb = list(B); op = rng.random()
            if op < 0.2 and len(nb) > 3: nb.pop(rng.randrange(len(nb)))
            elif op < 0.45: nb.append(rng.choice(nb))
            elif op < 0.6:
                v = rng.randint(1, 28); s = rng.choice([1, -1])
                if -s*v in nb: s = -s
                nb.append(s*v)
            elif op < 0.75:
                v = abs(rng.choice(nb)); nb = [(-x if abs(x) == v else x) for x in nb]
            else:
                i = rng.randrange(len(nb)); x = nb[i]; y = (1 if x > 0 else -1)*max(1, abs(x) + rng.choice([-1, 1]))
                if -y in nb: continue
                nb[i] = y
            if not pf(nb) or sum(abs(x) for x in nb) > WMAX: continue
            sc = score(nb)
            if sc is None: continue
            if sc[0] <= cur[0]: B, cur = nb, sc
            if cur[0] < 0: return (cur, sorted(B))
        if best is None or cur[0] < best[0][0]: best = (cur, sorted(B))
    return best
if __name__ == "__main__":
    with Pool(32) as pool: res = [r for r in pool.map(run, range(1000, 1032)) if r]
    res.sort(key=lambda r: r[0][0])
    for r in res[:4]: print(f"min normalized S-slack {r[0][0]:.4e} neg-p {r[0][1]} background {r[1]}")
    print("[S]", "VIOLATION FOUND" if res and res[0][0][0] < 0 else "no violation")
