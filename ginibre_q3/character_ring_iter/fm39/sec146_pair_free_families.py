from collections import defaultdict
import sys
def step(st, n, e):
    nx = defaultdict(int)
    for (s, t), v in st.items():
        for s2 in range(abs(s-n), s+n+1, 2): nx[s2, t] += v
        for t2 in range(abs(t-n), t+n+1, 2): nx[s, t2] += e*v
    return nx
def mstep(st, n):
    nx = defaultdict(int)
    for s, v in st.items():
        for s2 in range(abs(s-n), s+n+1, 2): nx[s2] += v
    return nx
for m in (5, 6, 7):
    for r in range(1, int(sys.argv[1])+1):
        st = {(0, 0): 1}; ms = {0: 1}; L = 0
        for n in range(3, m+1):
            for _ in range(r): st = step(st, n, -1); ms = mstep(ms, n); L += 1
        W = r*sum(range(3, m+1))
        best = None
        for p in range(m, W+1):
            if (W-p) % 2: continue
            d = (W-p)//2
            if d < 8 or m > d: continue
            if (L % 2 == 1) and p <= m: continue   # distinguished sign -: would pair with -p in background if p <= m
            g = st.get((p, 0), 0); mm = ms.get(p, 0)
            if mm and (best is None or g/mm < best[0]): best = (g/mm, p, g)
        if best: print(f"(-3..-{m})^{r}: L={L} W={W}  min ratio {best[0]:.4f} at p={best[1]}", flush=True)
