import random
from collections import defaultdict
exec(open(__file__.replace('sec147_mono_minus_new.py', 'sec147_mono_n.py')).read().split("random.seed")[0])
random.seed(13)
tot = viol = 0; worst = None; posfail = 0
for trial in range(4000):
    vals = random.sample(range(1, 13), random.randint(1, 5))
    lab = []
    for v in vals: lab += [-v]*random.randint(1, 4)          # all-minus background
    if random.random() < 0.5:                                  # optional single odd plus labels
        for v in random.sample([1,3,5,7,9,11], random.randint(1,2)):
            if -v not in lab: lab.append(v)
    n = random.choice([v for v in range(1, 15) if v not in vals and -v not in lab and v not in lab] or [99])
    if n == 99: continue
    g0 = gv(lab); g1 = step(step(g0, n, -1), n, -1)
    mx = max([abs(v) for v in lab] + [3]); W1 = sum(abs(v) for v in lab) + 2*n
    for p in range(mx, W1+1):
        if (W1 - p) % 2: continue
        if -p in lab or (p in lab): pass
        d = g1.get((p, 0), 0) - g0.get((p, 0), 0); tot += 1
        if d < 0:
            viol += 1
            if worst is None or d < worst[0]: worst = (d, sorted(lab), n, p)
print(f"all-minus backgrounds, add a NEW minus pair class (-n,-n): checked {tot}, violations {viol}, worst {worst}")
