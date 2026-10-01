# Exact two-core scan: background (+1)^a (-1)^e (+2)^b (eps m), distance >= 8, unsaturated (m <= delta), p >= max(m,6)
import sys
from collections import defaultdict
from fractions import Fraction
exec(open(__file__.replace('sec144_two_core_margin_scan.py', 'sec144_h1_margin_scan.py')).read().split("B = int")[0])
B = int(sys.argv[1]); NMAX = int(sys.argv[2]); MMAX = int(sys.argv[3])
rows = []
for m in range(3, MMAX+1):
  for sg in (1, -1):
    best = None; neg = 0
    for b in range(0, B+1):
      for a in range(0, NMAX+1):
        for e in range(0, NMAX+1):
            if a + e < 2: continue
            lab = [1]*a + [-1]*e + [2]*b + [sg*m]
            nminus = e + (1 if sg < 0 else 0)
            W = a + e + 2*b + m
            g = gvec(lab); mm = mvec(lab)
            for p in range(max(m, 6), W+1):
                if (W - p) % 2: continue
                delta = (W - p)//2
                if delta < 8 or m > delta: continue
                # distinguished sign fixed by parity: g_p is the FM3 value / 2 for that sign
                v = g.get((p, 0), 0); ref = mm.get(p, 0)
                if v < 0: neg += 1
                if ref:
                    r = Fraction(v, ref)
                    if best is None or r < best[0]: best = (r, a, e, b, p, delta)
    print(f"core {'+' if sg>0 else '-'}{m:2d}: negatives {neg}, min ratio {float(best[0]):.3e} at a={best[1]} e={best[2]} b={best[3]} p={best[4]} delta={best[5]}", flush=True)
