import random
from collections import defaultdict
exec(open(__file__.replace('sec147_mono_pair.py', 'sec147_mono_n.py')).read().split("random.seed")[0])
random.seed(12)
for (n, sgn) in ((1,-1),(3,-1),(4,-1),(5,-1),(7,-1),(1,1),(3,1),(5,1)):
    checked = viol = 0; worst = None
    for trial in range(1200):
        k = random.randint(1, 4)
        lab = [random.choice([1, -1])*random.randint(3, 12) for _ in range(k)]
        lab += [random.choice([1, -1])]*random.randint(0, 12) + [random.choice([2, -2])]*random.randint(0, 4)
        sg = {}; lab = [sg.setdefault(abs(v), 1 if v > 0 else -1)*abs(v) for v in lab]
        if sg.get(n, sgn) != sgn: continue
        g0 = gv(lab); g1 = step(step(g0, n, sgn), n, sgn)
        mx = max([abs(v) for v in lab] + [n, 3]); W1 = sum(abs(v) for v in lab) + 2*n
        for p in range(mx, W1+1):
            if (W1 - p) % 2: continue
            d = g1.get((p, 0), 0) - g0.get((p, 0), 0); checked += 1
            if d < 0:
                viol += 1
                if worst is None or d < worst[0]: worst = (d, sorted(lab), p)
    print(f"add ({'+' if sgn>0 else '-'}{n},{'+' if sgn>0 else '-'}{n}): checked {checked}, violations {viol}, worst {worst}")
