import argparse
import itertools
import time
from collections import defaultdict

parser = argparse.ArgumentParser(description="Exact FM-STR12b counterexample verifier")
parser.parse_args()

def stamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def cg(a, b):
    return range(abs(a - b), a + b + 1, 2)

def cg_table(word):
    d = {(0, 0): 1}
    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), v in d.items():
            for u in cg(r, n):
                q[u, s] += v
            for u in cg(s, n):
                q[r, u] += eps * v
        d = {k: v for k, v in q.items() if v}
    return d

def mul_poly(a, b):
    out = defaultdict(int)
    for (i, j), x in a.items():
        for (k, ell), y in b.items():
            out[i + k, j + ell] += x * y
    return {k: v for k, v in out.items() if v}

def U_polys(nmax):
    U = [[1], [0, 1]]
    for n in range(1, nmax):
        nxt = [0] * (n + 2)
        for k, v in enumerate(U[n]):
            nxt[k + 1] += v
        for k, v in enumerate(U[n - 1]):
            nxt[k] -= v
        U.append(nxt)
    return U[:nmax + 1]

def monomial_to_U(nmax):
    # x U_j = U_(j+1) + U_(j-1), with U_(-1)=0.
    out = [{0: 1}]
    for _ in range(nmax):
        nxt = defaultdict(int)
        for j, v in out[-1].items():
            nxt[j + 1] += v
            if j:
                nxt[j - 1] += v
        out.append({j: v for j, v in nxt.items() if v})
    return out

def polynomial_table(word):
    # Expand in x,y, then convert each monomial to U_r(x) U_s(y).
    nmax = sum(abs(z) for z in word)
    U = U_polys(nmax)
    M = monomial_to_U(nmax)
    poly = {(0, 0): 1}
    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        factor = {}
        for k, c in enumerate(U[n]):
            if c:
                factor[k, 0] = factor.get((k, 0), 0) + c
                factor[0, k] = factor.get((0, k), 0) + eps * c
        poly = mul_poly(poly, factor)
    table = defaultdict(int)
    for (i, j), c in poly.items():
        for r, a in M[i].items():
            for s, b in M[j].items():
                table[r, s] += c * a * b
    return {k: v for k, v in table.items() if v}

def layers(fa, fb, lam):
    d = defaultdict(int)
    for (r, s), a in fa.items():
        b = fb.get((r, s), 0)
        if b:
            d[lam[0] * r + lam[1] * s] += a * b
    return dict(sorted((k, v) for k, v in d.items() if v))

def prefixes(layer):
    out = {}
    total = 0
    for k, v in sorted(layer.items()):
        total += v
        out[k] = total
    return out

A = (-1, -1)
B = (1, 1, 2)
fa = cg_table(A)
fb = cg_table(B)
assert fa == polynomial_table(A)
assert fb == polynomial_table(B)
assert fa == {(0, 0): 2, (2, 0): 1, (1, 1): -2, (0, 2): 1}
assert fb == {
    (0, 0): 2, (2, 0): 3, (0, 2): 3, (4, 0): 1,
    (2, 2): 2, (1, 1): 4, (3, 1): 2, (1, 3): 2, (0, 4): 1
}
P = {k: v * fb.get(k, 0) for k, v in fa.items() if fb.get(k, 0)}
assert P == {(0, 0): 4, (2, 0): 3, (1, 1): -8, (0, 2): 3}

L10 = layers(fa, fb, (1, 0))
L12 = layers(fa, fb, (1, 2))
L1m1 = layers(fa, fb, (1, -1))
L11 = layers(fa, fb, (1, 1))
assert L10 == {0: 7, 1: -8, 2: 3}
assert L12 == {0: 4, 2: 3, 3: -8, 4: 3}
assert L1m1 == {-2: 3, 0: -4, 2: 3}
assert L11 == {0: 4, 2: -2}
assert prefixes(L10)[1] == -1
assert prefixes(L12)[3] == -1
assert prefixes(L1m1)[0] == -1
assert sum(P.values()) == 2

print(stamp(), "counterexample exact tables PASS", flush=True)
print("A coefficients", dict(sorted(fa.items())), flush=True)
print("B coefficients", dict(sorted(fb.items())), flush=True)
print("P=f_A*f_B", dict(sorted(P.items())), flush=True)
print("lambda=(1,0), layers", L10, "prefix T=1", prefixes(L10)[1], flush=True)
print("lambda=(1,2), layers", L12, "prefix T=3", prefixes(L12)[3], flush=True)
print("lambda=(1,-1), layers", L1m1, "prefix T=0", prefixes(L1m1)[0], flush=True)
print("lambda=(1,1), layers", L11, "all prefixes", prefixes(L11),
      "Phi", sum(P.values()), flush=True)

# Exact bounded scan: lengths 2..4, labels <=10, even minus count,
# proper splits counted once up to exchanging the two sides.
alphabet = tuple(range(-10, 0)) + tuple(range(1, 11))
checked = 0
for N in range(2, 5):
    found = None
    for word in itertools.combinations_with_replacement(alphabet, N):
        if sum(z < 0 for z in word) % 2:
            continue
        for mask in range(1, (1 << (N - 1))):
            left = tuple(word[i] for i in range(N) if mask >> i & 1)
            right = tuple(word[i] for i in range(N) if not (mask >> i & 1))
            checked += 1
            x = layers(cg_table(left), cg_table(right), (1, 0))
            pref = prefixes(x)
            if any(v < 0 for v in pref.values()):
                found = (word, left, right, pref)
                break
        if found:
            break
    print(stamp(), "bounded search length", N,
          "checked splits", checked, "first failure", found, flush=True)
    assert found is None
assert checked == 33605
print(stamp(), "bounded search through length 4 PASS", flush=True)
