# Exact H = 1 scan: background (+1)^a (-1)^e (+2)^b, all valid p; min ratio g_p / m_p(all) per b.
import sys
from collections import defaultdict
from fractions import Fraction
def gvec(labels):
    st = {(0, 0): 1}
    for x in labels:
        n = abs(x); e = 1 if x > 0 else -1; nx = defaultdict(int)
        for (s, t), v in st.items():
            for s2 in range(abs(s-n), s+n+1, 2): nx[s2, t] += v
            for t2 in range(abs(t-n), t+n+1, 2): nx[s, t2] += e*v
        st = nx
    return st
def mvec(labels):
    st = {0: 1}
    for x in labels:
        n = abs(x); nx = defaultdict(int)
        for s, v in st.items():
            for s2 in range(abs(s-n), s+n+1, 2): nx[s2] += v
        st = nx
    return st
B = int(sys.argv[1]); NMAX = int(sys.argv[2])
for b in range(0, B+1):
    best = None; neg = 0
    for a in range(0, NMAX+1):
        for e in range(0, NMAX+1):
            if (e % 2) or a + e < 2: continue       # even number of minus among background d's; distinguished sign by parity
            lab = [1]*a + [-1]*e + [2]*b
            W = a + e + 2*b
            g = gvec(lab); m = mvec(lab)
            for p in range(max(1, 2 if b else 1), W+1):
                if (W - p) % 2: continue
                delta = (W - p)//2
                if delta < 8 or p < 6: continue   # genuine core p >= 6 (labels <= 5 proved), distance >= 8
                v = g.get((p, 0), 0); mm = m.get(p, 0)
                if v < 0: neg += 1
                if mm:
                    r = Fraction(v, mm)
                    if best is None or r < best[0]: best = (r, a, e, p, delta)
    if best: print(f"b={b:2d}  negatives {neg}  min ratio {float(best[0]):.3e} at a={best[1]} e={best[2]} p={best[3]} delta={best[4]}", flush=True)
