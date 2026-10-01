# Test: g_p(a,e,b+1) >= g_p(a,e,b) for H = 1 backgrounds (+1)^a (-1)^e (+2)^b, all p >= 3.
from collections import defaultdict
import sys
exec(open(__file__.replace('sec145_h1_monotonicity_check.py', 'sec144_corner_scan.py')).read().split('P = int')[0])
NMAX = int(sys.argv[1]); BMAX = int(sys.argv[2])
viol = 0; checked = 0; worst = None
for a in range(0, NMAX+1):
    for e in range(0, NMAX+1):
        if a + e == 0: continue
        st = {(0, 0): 1}
        for _ in range(a): st = step(st, 1, 1)
        for _ in range(e): st = step(st, 1, -1)
        prev = {s: v for (s, t), v in st.items() if t == 0}
        for b in range(1, BMAX+1):
            st = step(st, 2, 1)
            cur = {s: v for (s, t), v in st.items() if t == 0}
            for p in range(3, a + e + 2*b + 1):
                if (a + e + 2*b - p) % 2: continue
                d = cur.get(p, 0) - prev.get(p, 0)
                checked += 1
                if d < 0:
                    viol += 1
                    if worst is None or d < worst[0]: worst = (d, a, e, b, p, prev.get(p,0), cur.get(p,0))
            prev = cur
print("checked", checked, "violations", viol, "worst", worst)
