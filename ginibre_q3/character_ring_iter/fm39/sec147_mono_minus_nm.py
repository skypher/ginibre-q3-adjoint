import random
from collections import defaultdict
exec(open(__file__.replace('sec147_mono_minus_nm.py', 'sec147_mono_n.py')).read().split("random.seed")[0])
random.seed(14)
tot = viol = 0; worst = None; ex = []
for trial in range(4000):
    vals = random.sample(range(1, 13), random.randint(1, 5))
    lab = []
    for v in vals: lab += [-v]*random.randint(1, 3)
    if random.random() < 0.5:
        for v in random.sample([1,3,5,7,9,11], random.randint(1,2)):
            if -v not in lab: lab.append(v)
    free = [v for v in range(1, 15) if v not in lab and -v not in lab]
    if len(free) < 2: continue
    n, m = random.sample(free, 2)
    g0 = gv(lab); g1 = step(step(g0, n, -1), m, -1)
    mx = max([abs(v) for v in lab] + [n, m, 3]); W1 = sum(abs(v) for v in lab) + n + m
    for p in range(mx, W1+1):
        if (W1 - p) % 2: continue
        d = g1.get((p, 0), 0) - g0.get((p, 0), 0); tot += 1
        if d < 0:
            viol += 1
            if len(ex) < 3: ex.append((sorted(lab), n, m, p, g0.get((p,0),0), g1.get((p,0),0)))
print(f"all-minus backgrounds, add two NEW distinct minus labels (-n,-m): checked {tot}, violations {viol}")
for e in ex: print("  example", e)
