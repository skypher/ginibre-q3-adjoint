"""FM-MECH106: FM3 for the whole cone follows from the pair reduction (FM-MECH102)
and two insertion-monotonicity lemmas on pair-free backgrounds (FM-SEC147):
 (M+) adding plus labels (an even label alone, or two odd labels) to a pair-free
      background B never decreases g_p(B), p >= max(labels, 3);
 (M-) for an all-minus background with at most one odd plus label, adding two
      minus labels never decreases g_p.
Descent: strip even plus labels one at a time and odd plus labels in pairs
(M+), then minus labels in pairs (M-), keeping the distinguished label p.  The
base lists have at most three factors: (p), (+-n, +-p) and (+o, -n, -p).
This script verifies that every base list has FM3 >= 0 (FM3 = 0, 2[n = p] and
2 m(o,n,p) respectively) for all labels <= 40, and replays the descent on 300
random pair-free lists, checking each step's inequality exactly."""
import random
from collections import defaultdict
def step(st, n, e):
    nx = defaultdict(int)
    for (s, t), v in st.items():
        for s2 in range(abs(s-n), s+n+1, 2): nx[s2, t] += v
        for t2 in range(abs(t-n), t+n+1, 2): nx[s, t2] += e*v
    return nx
def gp(lab, p):
    st = {(0, 0): 1}
    for x in lab: st = step(st, abs(x), 1 if x > 0 else -1)
    return st.get((p, 0), 0)
def fm3(word):
    st = {(0, 0): 1}
    for x in word: st = step(st, abs(x), 1 if x > 0 else -1)
    return st.get((0, 0), 0)
def m3(a, b, c):  # trivial multiplicity in U_a U_b U_c
    return int((a+b+c) % 2 == 0 and abs(a-b) <= c <= a+b)
cnt = 0
for p in range(1, 41):
    assert fm3([p]) == 0
    for n in range(1, 41):
        assert fm3([n, p]) == 2*(n == p) and fm3([-n, -p]) == 2*(n == p)
        cnt += 2
for o in range(1, 41, 2):
    for n in range(1, 41):
        for p in range(max(o, n), 41):
            assert fm3([o, -n, -p]) == 2*m3(o, n, p); cnt += 1
print("base lists: PASS,", cnt, "exact values")
random.seed(106); steps = 0
for trial in range(300):
    sg = {}
    lab = []
    for v in random.sample(range(1, 13), random.randint(2, 5)):
        sg[v] = random.choice([1, -1]); lab += [sg[v]*v]*random.randint(1, 3)
    W = sum(abs(v) for v in lab); mx = max(abs(v) for v in lab)
    ps = [p for p in range(max(mx, 3), W+1) if (W-p) % 2 == 0 and sg.get(p, 0) == 0]
    if not ps: continue
    p = random.choice(ps)
    cur = list(lab)
    # phase 1: strip plus labels
    while True:
        ev = [v for v in cur if v > 0 and v % 2 == 0]
        od = [v for v in cur if v > 0 and v % 2 == 1]
        if ev: rem = [ev[0]]
        elif len(od) >= 2: rem = od[:2]
        else: break
        smaller = list(cur)
        for r in rem: smaller.remove(r)
        assert gp(cur, p) >= gp(smaller, p), ("M+ step fails", cur, rem, p)
        cur = smaller; steps += 1
    # phase 2: strip minus labels in pairs
    while sum(1 for v in cur if v < 0) >= 2:
        mins = [v for v in cur if v < 0][:2]
        smaller = list(cur)
        for r in mins: smaller.remove(r)
        assert gp(cur, p) >= gp(smaller, p), ("M- step fails", cur, mins, p)
        cur = smaller; steps += 1
    assert len([v for v in cur if v < 0]) <= 1 and len([v for v in cur if v > 0]) <= 1
    assert gp(cur, p) >= 0
print("descent replay: PASS,", steps, "monotone steps on random pair-free lists")
