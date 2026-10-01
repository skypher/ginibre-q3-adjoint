# Independent check of the two Lemma R counterexamples (main agent, own two-spin DP).
from collections import defaultdict
def step(st, n, e):
    nx = defaultdict(int)
    for (a, b), v in st.items():
        for t in range(abs(a-n), a+n+1, 2): nx[t, b] += v
        for t in range(abs(b-n), b+n+1, 2): nx[a, t] += e*v
    return nx
def g(lab, p):
    st = {(0, 0): 1}
    for x in lab: st = step(st, abs(x), 1 if x > 0 else -1)
    return st.get((p, 0), 0)
B1 = [-1] + [-2]*2 + sum(([-n]*2 for n in range(3, 16, 2)), [])
B2 = [-1]*2 + sum(([n]*2 for n in range(3, 16, 2)), [])
for name, B, p, R in (("Minerva W=131", B1, 15, [-2, -2]), ("Juno W=128", B2, 16, [-1, -1])):
    sm = list(B)
    for r in R: sm.remove(r)
    a, b = g(B, p), g(sm, p)
    print(f"{name}: W={sum(abs(x) for x in B)} p={p} g_p(B)={a} g_p(B-R)={b} diff={a-b} -> {'VIOLATION' if a < b else 'ok'}")
