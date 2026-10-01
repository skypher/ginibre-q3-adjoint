from itertools import combinations
exec(open(__file__.replace('sec151_rule1_counterexample.py', 'sec151_rule1_adversarial.py')).read().split("def pf")[0])
B = [-1, -1, 2, 2, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11]; p = 11
W = sum(abs(v) for v in B); print("W", W, "parity ok", (W - p) % 2 == 0)
g = gv(B).get((p, 0), 0); print("g_p(B) =", g)
rems = [(v,) for v in sorted(set(B), key=abs) if abs(v) % 2 == 0] + sorted(set(tuple(sorted(c, key=abs)) for c in combinations(B, 2) if (abs(c[0]) + abs(c[1])) % 2 == 0), key=lambda t: [abs(x) for x in t])
ok = []
for R in rems:
    sm = list(B)
    for r in R: sm.remove(r)
    v = gv(sm).get((p, 0), 0)
    if g >= v: ok.append((R, v))
    if R == (-1, -1): print("rule removal (-1,-1): g_p(B - R) =", v, "->", "VIOLATION" if v > g else "ok")
print("monotone removals available:", [r for r, _ in ok][:20], "... total", len(ok))
