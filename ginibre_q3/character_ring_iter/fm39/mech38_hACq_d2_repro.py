"""FM-MECH38 (astra_max_ceres): H_AC_q proved for d <= 2 (uniform q-certificates); fusion-trajectory grouping obstructions."""
import argparse, ast
from pathlib import Path
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations
from math import comb

argparse.ArgumentParser(
    description="FM-MECH38 exact q-certificates and fusion-path census"
).parse_args()

def definitions(name):
    p = Path("ginibre_q3/character_ring_iter/fm39") / name
    tree = ast.parse(p.read_text())
    ns = {}
    nodes = [
        n for n in tree.body
        if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))
    ]
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), "exec"), ns)
    return ns

q = definitions("sec77_qgauss_repro.py")
add, mul, qb, fac, lin = (
    q[k] for k in ("add", "mul", "qbinom", "qfactorial", "linearization")
)
base = definitions("mech34_hAC_d2_repro.py")["certificate"]

def qi(n):
    return (1,)*n if n > 0 else (0,)

def scale(p, a):
    return tuple(a*x for x in p)

def psum(items):
    z = (0,)
    for p in items:
        z = add(z, p)
    return z

@lru_cache(None)
def row(mu):
    cur = {0: (1,)}
    for b in mu:
        nxt = {}
        for a, p in cur.items():
            for k in range(min(a, b)+1):
                j = a+b-2*k
                nxt[j] = add(nxt.get(j, (0,)), mul(p, lin(a, b, k)))
        cur = nxt
    return cur

def JK(mu):
    J = K = (0,)
    pref = []
    s = 0
    for x in mu:
        pref.append(s)
        J = add(J, mul(qi(s), qi(x)))
        K = add(K, mul(mul(qb(s, 2), qb(x, 2)), fac(2)))
        s += x
    for i, j in combinations(range(len(mu)), 2):
        term = (1,)
        for z in (pref[i], mu[i], pref[j]-2, mu[j]):
            term = mul(term, qi(z))
        K = add(K, term)
    return J, K

def square(p):
    out = {}
    for x, a in p.items():
        for y, b in p.items():
            out[x ^ y] = out.get(x ^ y, F(0)) + a*b
    return out

def certificate(mu, d):
    L = len(mu)
    t, v = mu.count(1), mu.count(2)
    I = [1 << i for i, x in enumerate(mu) if x == 1]
    Z = [1 << i for i, x in enumerate(mu) if x == 2]
    E = {x: 1 for x in I}
    A = {i ^ j: 1 for i, j in combinations(I, 2)}
    B = {x: 1 for x in Z}

    if d == 1:
        J = JK(mu)[0]
        rho = add(J, (-F(t, 2),))
        assert all(x >= 0 for x in rho)
        out = {S: (F(x, 2),) for S, x in square(E).items()}
        out[0] = add(out.get(0, (0,)), rho)
        return [out.get(S, (0,)) for S in range(1 << L)]

    removed = [i for i, x in enumerate(mu) if x == 1][:2]
    J = (
        JK(tuple(x for i, x in enumerate(mu) if i not in removed))[0]
        if t >= 2 else (0,)
    )
    K = JK(mu)[1]
    R = add(K, scale(J, -F(t, 2)))
    assert all(x >= 0 for x in R)
    K1 = K[1] if len(K) > 1 else 0
    assert K1 == 2*L*L-L*t-3*L-t-v+1
    J1 = J[1] if len(J) > 1 else 0
    if t >= 2:
        assert J1 == 2*L-t-4

    atoms = []
    def atom(c, p):
        assert c >= 0, (mu, c)
        if c and p:
            atoms.append((c, p))

    if t <= 1:
        atom(F(1, 2), B)
    elif t <= 3:
        atom(F(1, 2), A | B)
        atom(F(J1-t+2, 2), E)
    else:
        if v:
            for z in Z:
                atom(F(1, 6*v), A | {z: 3*v})
        else:
            atom(F(1, 6), A)
        atom(F(1, 2), B)
        atom((F(J1)-F(t-2, 3))/2, E)

    atom(
        F(K1)-sum(c*sum(x*x for x in p.values()) for c, p in atoms),
        {0: 1}
    )
    one = {}
    for c, p in atoms:
        for S, x in square(p).items():
            one[S] = one.get(S, F(0)) + c*x

    zero, _ = base(mu)
    sq = square(E)
    out = []
    for S in range(1 << L):
        high = add(
            scale(J, F(sq.get(S, 0), 2)),
            R if S == 0 else (0,)
        )
        out.append(add(
            (zero.get(S, 0), one.get(S, 0)),
            (0, 0)+high[2:]
        ))
    return out

counts = [0, 0]
entries = [0, 0]
for L in range(1, 7):
    for mu in combinations_with_replacement(range(1, 4), L):
        J, K = JK(mu)
        assert row(mu).get(sum(mu)-2, (0,)) == J
        assert row(mu).get(sum(mu)-4, (0,)) == K
        for d in (1, 2):
            n = sum(mu)-2*d
            if n < 1:
                continue
            target = certificate(mu, d)
            for S, p in enumerate(target):
                a = tuple(mu[i] for i in range(L) if S >> i & 1)
                b = tuple(mu[i] for i in range(L) if not (S >> i & 1))
                expected = mul(
                    row(a).get(0, (0,)), row(b).get(n, (0,))
                )
                assert add(p, scale(expected, -1)) == (0,)
                entries[d-1] += 1
            counts[d-1] += 1

print("d=1:", counts[0], "lists;", entries[0], "entries")
print("d=2:", counts[1], "lists;", entries[1], "entries")

# Positive product identities used in Lemma 3.
dchecks = 0
for M in [0]+list(range(2, 9)):
    for t in range(2, 10):
        if M == 0 and t < 5:
            continue
        T1 = psum(qi(M+i) for i in range(t-2))
        T2 = psum(
            mul(qi(M+i), qi(M+j-2))
            for i, j in combinations(range(t), 2)
        )
        D = add(T2, scale(T1, -F(t, 2)))
        if M == 0:
            r = t-3
            P = psum(
                mul(qi(k), psum(add(qi(i), (-1,))
                                for i in range(1, k+2)))
                for k in range(1, r+1)
            )
            P = add(P, scale(psum(
                add(qi(j), scale(qi(i), -1))
                for i, j in combinations(range(1, r+1), 2)
            ), F(1, 2)))
        else:
            P = psum(
                mul(qi(M+j-2), psum(add(qi(M+i), (-1,))
                                     for i in range(j)))
                for j in range(1, t)
            )
            P = add(P, qi(M-1))
            P = add(P, scale(T1, F(1, 2)))
            P = add(P, scale(psum(
                add(qi(M+j), scale(qi(M+i), -1))
                for i, j in combinations(range(t-2), 2)
            ), F(1, 2)))
        assert add(D, scale(P, -1)) == (0,)
        assert all(x >= 0 for x in P)
        dchecks += 1
print("positive tail identities:", dchecks)

def groups(labels):
    L = len(labels)
    N = 1 << (L-1)
    out = {}
    for S in range(N):
        colours = [(S >> i) & 1 for i in range(L)]
        rem = [
            sum(n for n, c in zip(labels, colours) if c == z)
            for z in (0, 1)
        ]
        if any(s % 2 for s in rem):
            continue
        cur = {(0, 0, (0,)): (1,)}
        for b, c in zip(labels, colours):
            rem[c] -= b
            nxt = {}
            for (u, v, path), p in cur.items():
                for k in range(min((u, v)[c], b)+1):
                    z = [u, v]
                    z[c] += b-2*k
                    if any(z[i] > rem[i] for i in (0, 1)):
                        continue
                    key = (z[0], z[1], path+(sum(z),))
                    nxt[key] = add(
                        nxt.get(key, (0,)),
                        mul(p, lin((u, v)[c], b, k))
                    )
            cur = nxt
        for (u, v, path), p in cur.items():
            assert u == v == 0
            if path not in out:
                out[path] = [(0,)]*N
            out[path][S] = add(out[path][S], p)
    return out

def walsh(a):
    a = list(a)
    h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            for j in range(i, i+h):
                a[j], a[j+h] = a[j]+a[j+h], a[j]-a[j+h]
        h *= 2
    return a

def test(tab):
    bad = None
    number = max(map(len, tab))
    for d in range(number):
        a = walsh([p[d] if d < len(p) else 0 for p in tab])
        if min(a) < 0 and bad is None:
            bad = (d, a.index(min(a)), min(a))
    return number, bad

raw = [0, 0, 0, 0]
rawfirst = None

def averaged(mu):
    global rawfirst
    L = len(mu)
    N = 1 << L
    distinct = sorted(set(mu))
    width = [mu.count(x) for x in distinct]
    orders = sorted(set(permutations(mu)))
    acc = {}

    for word in orders:
        raw[0] += 1
        G = groups(word)
        total = [(0,)]*(N//2)
        for path in sorted(G):
            tab = G[path]
            number, bad = test(tab)
            raw[1] += 1
            raw[2] += number
            raw[3] += bad is not None
            if bad is not None and rawfirst is None:
                rawfirst = (word, path)+bad
            for S, p in enumerate(tab):
                total[S] = add(total[S], p)
            if path not in acc:
                acc[path] = {}
            for S in range(N):
                key = tuple(
                    sum(1 for i, z in enumerate(word)
                        if z == x and S >> i & 1)
                    for x in distinct
                )
                p = tab[S if S < N//2 else (N-1) ^ S]
                acc[path][key] = add(acc[path].get(key, (0,)), p)

        for S, p in enumerate(total):
            a = tuple(word[i] for i in range(L) if S >> i & 1)
            b = tuple(word[i] for i in range(L) if not (S >> i & 1))
            assert p == mul(row(a).get(0, (0,)), row(b).get(0, (0,)))

    tables = {}
    for path, orbit in acc.items():
        tab = []
        for S in range(N//2):
            key = tuple(
                sum(1 for i, z in enumerate(mu)
                    if z == x and S >> i & 1)
                for x in distinct
            )
            den = len(orders)
            for w, k in zip(width, key):
                den *= comb(w, k)
            tab.append(tuple(F(x, den) for x in orbit.get(key, (0,))))
        tables[path] = tab
    return tables

avg = [0, 0, 0, 0]
first = witness = None
for L in range(1, 7):
    for mu in combinations_with_replacement(range(1, 4), L):
        if sum(mu) % 2:
            continue
        avg[0] += 1
        G = averaged(mu)
        total = [(0,)]*(1 << (L-1))
        for path in sorted(G):
            tab = G[path]
            number, bad = test(tab)
            avg[1] += 1
            avg[2] += number
            avg[3] += bad is not None
            if bad is not None and first is None:
                first = (mu, path)+bad
                witness = tab
            for S, p in enumerate(tab):
                total[S] = add(total[S], p)
        for S, p in enumerate(total):
            a = tuple(mu[i] for i in range(L) if S >> i & 1)
            b = tuple(mu[i] for i in range(L) if not (S >> i & 1))
            assert p == mul(row(a).get(0, (0,)), row(b).get(0, (0,)))
    print("through length", L, "averaged census", avg, flush=True)

assert raw == [546, 4757, 59266, 1501]
assert avg == [43, 3670, 45961, 343]
assert rawfirst == ((1,1,1,1), (0,1,2,1,0), 0, 3, -1)
assert first == (
    (1,1,1,1,2), (0,1,2,2,1,0), 0, 3, F(-1,15)
)

for S, p in enumerate(witness):
    wanted = (
        (F(1,5), F(3,5), F(3,5), F(1,5)) if S == 0 else
        ((F(2,15), F(2,15)) if S.bit_count() == 2 else (0,))
    )
    assert p == wanted
    repaired = add(p, (F(1,15),) if S == 0 else (0,))
    sq = 4 if S == 0 else (2 if S.bit_count() == 2 else 0)
    cert = add(
        scale((1,1), F(sq,15)),
        (0, F(1,3), F(3,5), F(1,5)) if S == 0 else (0,)
    )
    assert repaired == cert

print("raw census:", raw)
print("first raw failure:", rawfirst)
print("first averaged failure:", first)
print("all assertions passed")
