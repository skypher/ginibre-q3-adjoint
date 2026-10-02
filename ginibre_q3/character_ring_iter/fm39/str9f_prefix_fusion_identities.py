import sys
if any(arg in ("-h", "--help") for arg in sys.argv[1:]):
    print("FM-STR9f verifier: exact height-prefix fusion and split-move identities.")
    raise SystemExit(0)
from datetime import datetime, timezone
from itertools import product as cartesian_product
from collections import defaultdict

def log(*items):
    print(datetime.now(timezone.utc).isoformat(timespec="seconds"), *items, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def add_scaled(dst, src, scale=1):
    out = dict(dst)
    for key, value in src.items():
        out[key] = out.get(key, 0) + scale*value
    return {key: value for key, value in out.items() if value}

def multiply_factor(table, n, eps=1):
    out = defaultdict(int)
    for (r, s), value in table.items():
        for q in cg(r, n):
            out[q, s] += value
        for q in cg(s, n):
            out[r, q] += eps*value
    return {key: value for key, value in out.items() if value}

def character(word):
    row = {(0, 0): 1}
    for n, eps in word:
        row = multiply_factor(row, n, eps)
    return row

def F(n, eps):
    out = {}
    out[n, 0] = 1
    out[0, n] = out.get((0, n), 0) + eps
    return {key: value for key, value in out.items() if value}

def poly_product(P, Q):
    out = defaultdict(int)
    for (r, s), x in P.items():
        for (u, v), y in Q.items():
            for i in cg(r, u):
                for j in cg(s, v):
                    out[i, j] += x*y
    return {key: value for key, value in out.items() if value}

def prefix_pair(A, G, T):
    fa = character(A)
    return sum(value*G.get((r, s), 0)
               for (r, s), value in fa.items() if r+s <= T)

def degree(table):
    return max((r+s for r, s in table), default=0)

def pair_fusion_rhs(a, eps, b, eta):
    fused = {}
    for c in cg(a, b):
        fused = add_scaled(fused, F(c, eps*eta))
    flipped = poly_product(F(a, -eps), F(b, -eta))
    return add_scaled({k: 2*v for k, v in fused.items()}, flipped, -1)

def expansion3(labels, signs):
    a, b, c = labels
    e, f, g = signs
    out = {}
    for u in cg(a, b):
        for v in cg(u, c):
            out = add_scaled(out, F(v, e*f*g), 4)
    for u in cg(a, b):
        term = poly_product(F(u, -e*f), F(c, -g))
        out = add_scaled(out, term, -2)
    tail = poly_product(poly_product(F(a, -e), F(b, -f)), F(c, g))
    return add_scaled(out, tail, -1)

def expansion4(labels, signs):
    a, b, c, d = labels
    e, f, g, h = signs
    out = {}
    for u in cg(a, b):
        for v in cg(u, c):
            for w in cg(v, d):
                out = add_scaled(out, F(w, e*f*g*h), 8)
    for u in cg(a, b):
        for v in cg(u, c):
            out = add_scaled(out, poly_product(F(v, -e*f*g), F(d, -h)), -4)
    for u in cg(a, b):
        term = poly_product(poly_product(F(u, -e*f), F(c, -g)), F(d, h))
        out = add_scaled(out, term, -2)
    tail = poly_product(poly_product(poly_product(F(a, -e), F(b, -f)), F(c, g)), F(d, h))
    return add_scaled(out, tail, -1)

def move_formula(A, B, n, eps, T):
    fa, fb = character(A), character(B)
    value = 0
    for (q, s), x in fa.items():
        for r in cg(q, n):
            value += x*fb.get((r, s), 0)*(
                int(r+s <= T)-int(q+s <= T))
    for (r, q), x in fa.items():
        for s in cg(q, n):
            value += eps*x*fb.get((r, s), 0)*(
                int(r+s <= T)-int(r+q <= T))
    return value

def prefix_direct(A, B, T):
    return prefix_pair(A, character(B), T)

def main():
    log("start exact height-prefix identities")
    pair_checks = 0
    for a in range(0, 5):
        for b in range(0, 5):
            for e, f in cartesian_product((-1, 1), repeat=2):
                lhs = poly_product(F(a, e), F(b, f))
                rhs = pair_fusion_rhs(a, e, b, f)
                assert lhs == rhs, (a, b, e, f, lhs, rhs)
                pair_checks += 1
    log("pair fusion identity PASS", pair_checks, "exact polynomial checks")

    Atests = [[], [(1, 1)], [(1, -1)], [(1, -1), (2, -1)]]
    checks3 = checks4 = 0
    for labels in cartesian_product(range(1, 4), repeat=3):
        for signs in cartesian_product((-1, 1), repeat=3):
            direct = character(list(zip(labels, signs)))
            expanded = expansion3(labels, signs)
            assert direct == expanded, ("3-factor identity", labels, signs)
            for A in Atests:
                maxT = degree(poly_product(character(A), direct))
                for T in range(maxT+1):
                    lhs = prefix_direct(A, list(zip(labels, signs)), T)
                    rhs = prefix_pair(A, expanded, T)
                    assert lhs == rhs
                    checks3 += 1
        if labels[1:] == (1, 1):
            log("3-factor expansions", "first label", labels[0],
                "prefix pairings", checks3)
    log("3-factor expansion PASS", checks3, "prefix pairings")

    for labels in cartesian_product(range(1, 3), repeat=4):
        for signs in cartesian_product((-1, 1), repeat=4):
            direct = character(list(zip(labels, signs)))
            expanded = expansion4(labels, signs)
            assert direct == expanded, ("4-factor identity", labels, signs)
            for A in Atests:
                maxT = degree(poly_product(character(A), direct))
                for T in range(maxT+1):
                    lhs = prefix_direct(A, list(zip(labels, signs)), T)
                    rhs = prefix_pair(A, expanded, T)
                    assert lhs == rhs
                    checks4 += 1
        if labels[1:] == (1, 1, 1):
            log("4-factor expansions", "first label", labels[0],
                "prefix pairings", checks4)
    log("4-factor expansion PASS", checks4, "prefix pairings")

    words = [[], [(1, 1)], [(1, -1)], [(2, 1)], [(2, -1)],
             [(1, 1), (2, -1)], [(1, -1), (2, -1)]]
    move_checks = 0
    for A in words:
        for B in words:
            for n in range(1, 4):
                for eps in (-1, 1):
                    maxT = degree(character(A+[(n, eps)]+B)) + 1
                    for T in range(maxT+1):
                        direct = (prefix_direct(A+[(n, eps)], B, T)
                                  - prefix_direct(A, B+[(n, eps)], T))
                        formula = move_formula(A, B, n, eps, T)
                        assert direct == formula, (A, B, n, eps, T, direct, formula)
                        move_checks += 1
    log("split-move commutator PASS", move_checks, "exact checks")

    A = [(1, -1)]*3
    B = [(2, -1), (2, -1), (1, -1)]
    T = 3
    Q = sum(prefix_direct(A, [(c, 1), (1, -1)], T) for c in (0, 2, 4))
    parent = prefix_direct(A, B, T)
    flipped = prefix_direct(A, [(2, 1), (2, 1), (1, -1)], T)
    assert (parent, Q, flipped, parent-Q) == (28, 40, 52, -12)
    assert parent == 2*Q-flipped
    all_signs = [eps for _, eps in A+B]
    assert len(A) == len(B) == 3 and sum(eps < 0 for eps in all_signs) % 2 == 0
    assert all(not (n, -eps) in set(A+B) for n, eps in set(A+B))
    log("pair-free signed-remainder witness", "A=(-1,-1,-1)",
        "B=(-2,-2,-1)", "T=3", "parent=28", "fused-child-sum=40",
        "flipped=52", "remainder=-12")

    dplus = (prefix_direct([(2, 1), (1, -1)], [(1, -1)], 1)
             - prefix_direct([(2, 1)], [(1, -1), (1, -1)], 1))
    dminus = (prefix_direct([(1, -1), (1, -1)], [(2, 1)], 1)
              - prefix_direct([(1, -1)], [(2, 1), (1, -1)], 1))
    assert (dplus, dminus) == (2, -2)
    log("split-move sign witness", "same labels and T=1", "differences=+2,-2")
    log("PASS FM-STR9f verifier")

main()
