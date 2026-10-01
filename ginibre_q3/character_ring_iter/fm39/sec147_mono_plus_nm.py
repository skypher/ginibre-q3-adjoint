import random
from collections import defaultdict
exec(open(__file__.replace('sec147_mono_plus_nm.py', 'sec147_mono_n.py')).read().split("random.seed")[0])
random.seed(15)
tot = viol = 0; ex = []
for trial in range(5000):
    k = random.randint(0, 4)
    lab = [random.choice([1, -1])*random.randint(3, 12) for _ in range(k)]
    lab += [random.choice([1, -1])]*random.randint(0, 12) + [random.choice([2, -2])]*random.randint(0, 4)
    sg = {}; lab = [sg.setdefault(abs(v), 1 if v > 0 else -1)*abs(v) for v in lab]
    if not lab: continue
    cand = [v for v in range(1, 16, 2) if sg.get(v, 1) > 0]
    if len(cand) < 2: continue
    n, m = random.sample(cand, 2)
    g0 = gv(lab); g1 = step(step(g0, n, 1), m, 1)
    mx = max([abs(v) for v in lab] + [n, m, 3]); W1 = sum(abs(v) for v in lab) + n + m
    for p in range(mx, W1+1):
        if (W1 - p) % 2: continue
        d = g1.get((p, 0), 0) - g0.get((p, 0), 0); tot += 1
        if d < 0:
            viol += 1
            if len(ex) < 3: ex.append((sorted(lab), n, m, p, g0.get((p,0),0), g1.get((p,0),0)))
print(f"pair-free backgrounds, add two distinct odd plus labels (+n,+m): checked {tot}, violations {viol}")
for e in ex: print("  example", e)
