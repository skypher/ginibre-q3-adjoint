"""FM-MECH102: the pair reduction.  For every list containing +n and -n,
FM3(rest, +n, -n) = sum_(k=1..n) FM3(rest, -2k), because
(U_n(x) + U_n(y))(U_n(x) - U_n(y)) = U_n(x)^2 - U_n(y)^2
   = sum_(k=1..n) (U_(2k)(x) - U_(2k)(y)).
Each summand has one factor fewer and the same minus parity, so by induction on
the number of factors FM3 for all lists follows from FM3 for pair-free lists.
This script checks the polynomial identity symbolically for n <= 12 and the
FM3 identity on 2,000 random lists with an independent fusion evaluator."""
import random
from collections import defaultdict
import sympy as sp

x, y = sp.symbols("x y")
def U(n, t):
    a, b = sp.Integer(1), t
    if n == 0: return a
    for _ in range(n-1): a, b = b, sp.expand(t*b - a)
    return b
for n in range(1, 13):
    lhs = sp.expand((U(n, x) + U(n, y))*(U(n, x) - U(n, y)))
    rhs = sp.expand(sum(U(2*k, x) - U(2*k, y) for k in range(1, n+1)))
    assert lhs == rhs
print("polynomial identity: PASS for n <= 12")

def fm3(word):
    st = {(0, 0): 1}
    for sn in word:
        n = abs(sn); s = 1 if sn > 0 else -1; nx = defaultdict(int)
        for (a, b), v in st.items():
            for t in range(abs(a-n), a+n+1, 2): nx[t, b] += v
            for t in range(abs(b-n), b+n+1, 2): nx[a, t] += s*v
        st = nx
    return st.get((0, 0), 0)
random.seed(102); ok = 0
for _ in range(2000):
    n = random.randint(1, 8)
    rest = [random.choice([1, -1])*random.randint(1, 9) for _ in range(random.randint(1, 7))]
    if sum(1 for v in rest if v < 0) % 2 == 0: rest.append(-random.randint(1, 3))
    # rest has odd minus count, so rest + (+n, -n) has even minus count
    lhs = fm3(rest + [n, -n]); rhs = sum(fm3(rest + [-2*k]) for k in range(1, n+1))
    assert lhs == rhs, (rest, n, lhs, rhs); ok += 1
print("FM3 pair identity: PASS on", ok, "random lists")
