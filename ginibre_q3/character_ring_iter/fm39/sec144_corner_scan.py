from collections import defaultdict
import math, sys
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
P = int(sys.argv[1])
for a in range(4, int(sys.argv[2])+1, 4):
    st = {(0, 0): 1}; ms = {0: 1}
    for _ in range(a): st = step(st, 1, 1); ms = mstep(ms, 1)
    for _ in range(a): st = step(st, 1, -1); ms = mstep(ms, 1)
    best = (1e9, -1); vals = []
    for b in range(0, 2*a + 12):
        if b: st = step(st, 2, 1); ms = mstep(ms, 2)
        W = 2*a + 2*b
        if (W - P) % 2 == 0 and (W - P)//2 >= 8:
            g = st.get((P, 0), 0); m = ms.get(P, 0)
            r = g / m if m else float('inf')
            if g < 0: print("NEGATIVE", a, b, P, g)
            if r < best[0]: best = (r, b)
    print(f"a=e={a:3d}: min over b of g_{P}/m_{P}(all) = {best[0]:.3e} at b={best[1]}  (log10 {math.log10(best[0]):.2f})", flush=True)
