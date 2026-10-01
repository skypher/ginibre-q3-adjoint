# Test: adding an even plus label +n (n = 4, 6, 8) never decreases g_p, for p >= max(labels, n, 3), pair-free kept.
import random
from collections import defaultdict
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
random.seed(11)
for n in (4, 6, 8):
    checked = viol = 0; worst = None
    for trial in range(1500):
        k = random.randint(1, 4)
        lab = [random.choice([1, -1])*random.randint(3, 12) for _ in range(k)]
        lab += [random.choice([1, -1])]*random.randint(0, 12) + [random.choice([2, -2])]*random.randint(0, 4)
        sg = {}; lab = [sg.setdefault(abs(v), 1 if v > 0 else -1)*abs(v) for v in lab]   # pair-free
        if sg.get(n, 1) < 0: continue           # adding +n must keep pair-free
        g0 = gv(lab); g1 = step(g0, n, 1)
        mx = max([abs(v) for v in lab] + [n, 3]); W1 = sum(abs(v) for v in lab) + n
        for p in range(mx, W1+1):
            if (W1 - p) % 2: continue
            d = g1.get((p, 0), 0) - g0.get((p, 0), 0); checked += 1
            if d < 0:
                viol += 1
                if worst is None or d < worst[0]: worst = (d, sorted(lab), p)
    print(f"+{n}: checked {checked}, violations {viol}, worst {worst}")
