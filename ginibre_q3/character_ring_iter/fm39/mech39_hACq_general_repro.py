"""FM-MECH39 (astra_max_ceres): general contraction formula; H_AC_q for mu_i >= d (all d); d = 3 with at most one label 1; origin-only transfer obstruction; (1^4,2^2) repair.  Run with --part proofs and --part transfers."""
import argparse, ast
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from collections import Counter
from math import comb, lcm

parser = argparse.ArgumentParser(description="FM-MECH39 exact verifier")
parser.add_argument("--part", choices=("proofs", "transfers"), default="proofs")
args = parser.parse_args()
root = Path("ginibre_q3/character_ring_iter/fm39")

def definitions(name):
    path = root / name
    tree = ast.parse(path.read_text())
    nodes = [x for x in tree.body
             if isinstance(x, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    env = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), env)
    return env

old = definitions("mech38_hACq_d2_repro.py")
q = definitions("sec77_qgauss_repro.py")
old.update(add=q["add"], mul=q["mul"], lin=q["linearization"],
           qb=q["qbinom"], fac=q["qfactorial"],
           base=definitions("mech34_hAC_d2_repro.py")["certificate"])
add, mul, lin, qb, fac = [old[k] for k in ("add", "mul", "lin", "qb", "fac")]
row, scale, square, qi, psum = [
    old[k] for k in ("row", "scale", "square", "qi", "psum")]
b, beta, gamma = (1, 1), (1, 3, 3, 1), (1, 2, 2, 1)

def difference(p, q):
    return add(p, scale(q, -1))

def target(mu, n):
    out = []
    for S in range(1 << len(mu)):
        left = tuple(x for i, x in enumerate(mu) if S >> i & 1)
        right = tuple(x for i, x in enumerate(mu) if not (S >> i & 1))
        out.append(mul(row(left).get(0, (0,)), row(right).get(n, (0,))))
    return out

def top(mu, d):
    prefix = [0]
    for x in mu:
        prefix.append(prefix[-1] + x)
    @lru_cache(None)
    def visit(i, used):
        if i == len(mu):
            return (int(used == d),)
        h = prefix[i] - 2*used
        return psum(mul(lin(h, mu[i], k), visit(i+1, used+k))
                    for k in range(min(h, mu[i], d-used)+1))
    return visit(0, 0)

def partial_matchings(mu, d):
    owner = tuple(i for i, x in enumerate(mu) for _ in range(x))
    count = Counter()
    def visit(todo, left, pairs, opens):
        if left == 0:
            opens += todo
            crossing = sum(a < c < z < w or c < a < w < z
                           for (a, z), (c, w) in combinations(pairs, 2))
            covered = sum(a < x < z for a, z in pairs for x in opens)
            count[crossing+covered] += 1
            return
        if len(todo) < 2*left:
            return
        a, rest = todo[0], todo[1:]
        visit(rest, left, pairs, opens+(a,))
        for z in rest:
            if owner[a] != owner[z]:
                visit(tuple(x for x in rest if x != z), left-1,
                      pairs+((a, z),), opens)
    visit(tuple(range(len(owner))), d, (), ())
    return tuple(count[k] for k in range(max(count, default=0)+1))

def expand(mu, origin, atoms):
    cost = psum(scale(c, sum(x*x for x in p.values())) for c, p in atoms)
    residual = difference(origin, cost)
    assert all(x >= 0 for x in residual), (mu, residual)
    out = [(0,)] * (1 << len(mu))
    for c, p in atoms + [(residual, {0: 1})]:
        assert all(x >= 0 for x in c), (mu, c)
        for S, x in square(p).items():
            out[S] = add(out[S], scale(c, x))
    return out

def large_minimum(mu, d):
    E = {1 << i: 1 for i, x in enumerate(mu) if x == d}
    return expand(mu, top(mu, d), [(scale(fac(d), F(1, 2)), E)])

def distance_three(mu):
    t, v, u = mu.count(1), mu.count(2), mu.count(3)
    assert t <= 1
    n = sum(mu)-6
    if mu == (1, 2, 2, 2):
        previous = old["certificate"]((1, 1, 2, 2), 2)
        out = []
        for S in range(16):
            e = S >> 3 & 1
            T = ((S & 1) ^ e) | (e << 1)
            T |= (((S >> 1 & 1) ^ e) << 2)
            T |= (((S >> 2 & 1) ^ e) << 3)
            out.append(mul(b, previous[T]))
        return out
    B = {1 << i: 1 for i, x in enumerate(mu) if x == 2}
    Z = {1 << i: 1 for i, x in enumerate(mu) if x == 3}
    e = next((1 << i for i, x in enumerate(mu) if x == 1), 0)
    A = (0,)
    if v >= 2:
        remaining = list(mu)
        remaining.remove(2)
        remaining.remove(2)
        A = mul(b, top(tuple(remaining), 1))
    atoms = []
    if not t:
        atoms += [(scale(A, F(1, 2)), B), (scale(gamma, F(1, 2)), Z)]
    elif v == 0:
        atoms += [(scale(gamma, F(1, 2)), Z)]
    elif v == 1:
        if u:
            atoms += [(scale(gamma, F(1, 2)), B | {e ^ z: 1 for z in Z})]
    else:
        atoms += [(scale(gamma, F(1, 2)), B | {e ^ z: 1 for z in Z}),
                  (scale(difference(A, gamma), F(1, 2)), B)]
    for i, j, k in combinations(B, 3):
        atoms.append((scale(beta, F(1, 2)), {i: 1, j ^ k: 1}))
    return expand(mu, top(mu, 3), atoms)

def proofs():
    paths = 0
    for L in range(1, 5):
        for mu in combinations_with_replacement(range(1, 4), L):
            for d in range(1, 5):
                if sum(mu) >= 2*d:
                    assert partial_matchings(mu, d) == row(mu).get(sum(mu)-2*d, (0,))
                    paths += 1
    print("partial-matching identities:", paths)
    contractions = 0
    for L in range(1, 8):
        for mu in combinations_with_replacement(range(1, 4), L):
            for d in (3, 4):
                if sum(mu)-2*d >= 1:
                    assert top(mu, d) == row(mu).get(sum(mu)-2*d, (0,))
                    contractions += 1
    print("general contraction formula, d=3,4:", contractions)
    counts = entries = 0
    for d in range(1, 7):
        for L in range(1, 6):
            for mu in combinations_with_replacement(range(d, d+3), L):
                n = sum(mu)-2*d
                if n < 1:
                    continue
                assert large_minimum(mu, d) == target(mu, n)
                counts += 1
                entries += 1 << L
    print("minimum-label theorem:", counts, "lists;", entries, "entries")
    counts, entries = [0, 0], [0, 0]
    for L in range(1, 9):
        for t in (0, 1):
            for rest in combinations_with_replacement(range(2, 6), L-t):
                mu = (1,)*t + rest
                if sum(mu)-6 < 1:
                    continue
                assert distance_three(mu) == target(mu, sum(mu)-6)
                counts[t] += 1
                entries[t] += 1 << L
    print("d=3, counts of ones 0/1:", counts, "entries:", entries)
    def J(v):
        return top((2,)*v, 1)
    def K(v):
        return top((2,)*v, 2)
    def R(v):
        return difference(top((2,)*v, 3),
                          add(scale(mul(b, J(v-2)), F(v, 2)),
                              scale(beta, comb(v, 3))))
    assert R(3) == (0,)
    for v in range(3, 13):
        A = F(comb(v, 2))-F(v+1, 2)
        first = mul(mul(b, qi(2*v-4)), difference(K(v), scale(b, comb(v, 2))))
        bracket = scale(difference(qi(2*v-4), b), A+1)
        bracket = add(bracket, difference(qi(2*v-2), b))
        bracket = add(bracket, scale(psum(difference(qi(2*i), b)
                                          for i in range(1, v-2)), F(1, 2)))
        bracket = add(bracket, mul(difference(qb(2*v-2, 2), (1,)),
                                   psum(qi(2*i) for i in range(1, v))))
        rhs = add(first, mul(mul(b, b), bracket))
        assert difference(R(v+1), R(v)) == rhs
        assert all(x >= 0 for x in rhs)
    print("positive all-2 recurrence: 10 identities; PASS")

def transfers():
    old.update(raw=[0, 0, 0, 0], rawfirst=None, test=lambda tab: (0, None))
    counts, failures = [0, 0, 0], []
    for L in range(2, 8):
        for lam in combinations_with_replacement(range(1, 4), L):
            if sum(lam) % 2:
                continue
            d = (sum(lam)-2*max(lam))//2
            if d not in (3, 4):
                continue
            groups = old["averaged"](lam)
            degree = max(len(p) for tab in groups.values() for p in tab)
            available, needed = [F(0)]*degree, [F(0)]*degree
            for tab in groups.values():
                den = 1
                for p in tab:
                    for z in p:
                        den = lcm(den, z.denominator)
                integer = [[int(z*den) for z in p] for p in tab]
                for k in range(degree):
                    col = [p[k] if k < len(p) else 0 for p in integer]
                    available[k] += F(col[0], den)
                    col[0] = 0
                    needed[k] += F(-min(old["walsh"](col)), den)
            bad = [k for k in range(degree) if needed[k] > available[k]]
            counts[0] += 1
            counts[1] += degree
            counts[2] += len(bad)
            failures += [(lam, d, k, available[k], needed[k]) for k in bad]
        print("transfer census through L", L, counts, flush=True)
    assert counts == [20, 347, 4]
    assert failures == [
        ((1,1,1,1,2,2,2),3,0,F(13),F(4321,315)),
        ((1,1,1,1,2,3,3),3,0,F(14),F(647,45)),
        ((1,1,1,2,2,2,3),3,0,F(17),F(597,35)),
        ((1,1,2,2,2,3,3),4,0,F(25),F(8944,315))]
    print("origin-only failures:", failures)
    print("PASS")

def repaired_witness():
    mu = (1, 1, 1, 1, 2, 2)
    E = {1, 2, 4, 8}
    B = {16, 32}
    A = {i ^ j for i, j in combinations(E, 2)}
    H = {i ^ j ^ k for i, j, k in combinations(E, 3)}
    J = {i ^ j for i in E for j in B}
    Z = {48}
    data = [
        (F(1,4), F(3,2), F(1,8), F(2), F(5,24),
         F(12,5), F(13,48), F(52,15)),
        (F(5,12), F(3,2), F(27,40), F(10,9), F(19,120),
         F(140,19), F(17,16), F(764,57)),
        (F(1,8), F(3), F(9,16), F(4,3), F(5,48),
         F(12), F(55,24), F(88,3)),
        (F(1,48), F(6), F(3,20), F(5,3), F(1,60),
         F(55,2), F(19,8), F(2203,48))]
    def combine(P, Q, r):
        return {x: F(x in P) + r*F(x in Q) for x in P | Q}
    out = [(0,)]*64
    for k, (y, r, a, s, c, t, e, rho) in enumerate(data):
        atoms = [(y, combine(J,H,r)), (a,combine(A,B,s)),
                 (c,combine(A,Z,t)), (e,{x:1 for x in E}), (rho,{0:1})]
        for w, p in atoms:
            assert w >= 0 and all(v >= 0 for v in p.values())
            for S, v in square(p).items():
                out[S] = add(out[S], (0,)*k+(w*v,))
    for k, w, rho in [(4,F(3,2),56),(5,F(1,2),41)]:
        for S, v in square({x:1 for x in E}).items():
            out[S] = add(out[S], (0,)*k+(w*v,))
        out[0] = add(out[0], (0,)*k+(rho,))
    out[0] = add(out[0], (0,)*6+(24,11,4,1))
    assert out == target(mu, 2)
    print("first failed transfer list: q-positive repair, 64 polynomial entries; PASS")

if args.part == "proofs":
    proofs()
    repaired_witness()
else:
    transfers()
