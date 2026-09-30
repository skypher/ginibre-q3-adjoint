"""FM-MECH34 (astra_max_ceres): H_AC for every list with a distinguished label n = sum(mu) - 4;
explicit nonnegative autocorrelation factors (Theorem 1), checked against direct fusion tables."""
import argparse, random
from math import comb
from fractions import Fraction as Q
from itertools import combinations

parser = argparse.ArgumentParser()
parser.parse_args()

def add_scaled(dst, src, c):
    for x, v in src.items():
        dst[x] = dst.get(x, Q(0)) + c*v

def square(p):
    out = {}
    for x, a in p.items():
        for y, b in p.items():
            out[x ^ y] = out.get(x ^ y, Q(0)) + a*b
    return out

def fusion_rows(labels):
    rows = [{0: 1}]
    for n in labels:
        old = rows[:]
        for row in old:
            out = {}
            for j, v in row.items():
                for q in range(abs(j-n), j+n+1, 2):
                    out[q] = out.get(q, 0) + v
            rows.append(out)
    return rows

def certificate(mu):
    L = len(mu)
    t, v = mu.count(1), mu.count(2)
    w = L-t-v
    I = [1 << i for i, n in enumerate(mu) if n == 1]
    J = [1 << i for i, n in enumerate(mu) if n == 2]
    E = {x: Q(1) for x in I}
    A = {x ^ y: Q(1) for x, y in combinations(I, 2)}
    B = {x: Q(1) for x in J}
    atoms = []

    def atom(c, p):
        assert c >= 0 and all(z >= 0 for z in p.values())
        if c and p:
            atoms.append((c, p))

    if t <= 1:
        atom(Q(1, 2), B)
        rho = Q(comb(L, 2)-t) - Q(v, 2)

    elif t <= 3:
        p = A.copy()
        add_scaled(p, B, Q(1))
        atom(Q(1, 2), p)
        gamma = v+w-1
        atom(Q(gamma, 2), E)
        rho = (Q(comb(L, 2)-t)
               - Q(comb(t, 2)+v, 2) - Q(gamma*t, 2))

    else:
        if v == 0:
            atom(Q(1, 3), A)
            z = Q(0)
        elif v <= 3:
            for x in J:
                p = A.copy()
                p[x] = Q(3*v, 2)
                atom(Q(1, 3*v), p)
            atom(Q(1, 2), B)
            z = Q(3*v*v, 4) + Q(v, 2)
        else:
            for x, y in combinations(J, 2):
                p = A.copy()
                p[x] = p[y] = Q(v-1)
                atom(Q(v, 4*(v-1)*comb(v, 2)), p)
            atom(Q(v-4, 12*(v-1)), A)
            z = Q(v*(v-1), 2)

        gamma = Q(t-5, 3)+v+w
        atom(gamma/2, E)
        rho = (Q(comb(L, 2)-t)
               - Q(comb(t, 2), 3) - z - gamma*t/2)

    atom(rho, {0: Q(1)})
    out = {}
    for c, p in atoms:
        add_scaled(out, square(p), c)
    return out, rho

cases = []
for L in range(1, 10):
    for t in range(L+1):
        for v in range(L-t+1):
            mu = (1,)*t + (2,)*v + (3,)*(L-t-v)
            if sum(mu) >= 5:
                cases.append(mu)

rng = random.Random(3401)
for _ in range(200):
    mu = tuple(rng.randrange(1, 9)
               for _ in range(rng.randrange(1, 11)))
    if sum(mu) >= 5:
        cases.append(mu)

entries = 0
least = None
for mu in cases:
    n = sum(mu)-4
    cert, rho = certificate(mu)
    rows = fusion_rows(mu)
    full = len(rows)-1
    for x in range(full+1):
        value = rows[x].get(0, 0)*rows[full ^ x].get(n, 0)
        assert cert.get(x, 0) == value
        entries += 1
    least = rho if least is None else min(least, rho)

print(len(cases), "lists;", entries, "direct fusion entries; PASS")
print("least residual coefficient:", least)
