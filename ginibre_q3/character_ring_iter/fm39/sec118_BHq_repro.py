"""FM-SEC118 (luna_max_venus): BH_q census (Bernstein-in-s tables of the braided family are H_AC_q); stages: smallcert, sixcert A B, boundary7, boundary9."""
from itertools import combinations, product
from fractions import Fraction as Q
from math import comb
from scipy.optimize import linprog
from sympy import Matrix, Rational
import sys

def add(a, b, scale=1):
    z = a.copy()
    for e, c in b.items():
        z[e] = z.get(e, 0) + scale*c
    return {e: c for e, c in z.items() if c}

def shift(a, dq=0, ds=0):
    return {(i+dq, j+ds): c for (i, j), c in a.items()}

def qmul(a, n):
    z = {}
    for k in range(n):
        z = add(z, shift(a, dq=k))
    return z

# Full Walsh transform of f_(q,s), as a polynomial in q and s.
def fhat(labels, T):
    def ann(col, V):
        out = {}
        for word, poly in V.items():
            factor = poly
            for pos, other in enumerate(word):
                if other == col:
                    w = word[:pos] + word[pos+1:]
                    out[w] = add(out.get(w, {}), factor)
                factor = shift(factor, dq=(other == col), ds=(other != col))
                if not factor:
                    break
        return {w: p for w, p in out.items() if p}

    def X(col, V, maxdeg):
        out = ann(col, V)
        for word, poly in V.items():
            w = (col,) + word
            if len(w) <= maxdeg:
                out[w] = add(out.get(w, {}), poly)
        return {w: p for w, p in out.items() if p}

    def H(n, col, V, maxdeg):
        prev = {}
        cur = {w: p.copy() for w, p in V.items()}
        for k in range(n):
            nxt = X(col, cur, maxdeg)
            if k:
                for w, p in prev.items():
                    nxt[w] = add(nxt.get(w, {}), qmul(p, k), -1)
            prev, cur = cur, {w: p for w, p in nxt.items() if p}
        return cur

    V = {(): {(0, 0): 1}}
    for i in range(len(labels)-1, -1, -1):
        maxdeg = sum(labels[:i])
        sign = -1 if (T >> i) & 1 else 1
        A = H(labels[i], 0, V, maxdeg + labels[i])
        B = H(labels[i], 1, V, maxdeg + labels[i])
        nxt = {}
        for w, p in A.items():
            if len(w) <= maxdeg:
                nxt[w] = add(nxt.get(w, {}), p)
        for w, p in B.items():
            if len(w) <= maxdeg:
                nxt[w] = add(nxt.get(w, {}), p, sign)
        V = {w: p for w, p in nxt.items() if p}
    return V.get((), {})

def wht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2*h):
            for j in range(i, i+h):
                x, y = v[j], v[j+h]
                v[j], v[j+h] = x+y, x-y
        h *= 2
    return v

def bernstein_data(labels):
    L = len(labels)
    G = 1 << (L-1)
    # The dual group consists of even-cardinality masks.
    raw = [
        fhat(labels, x | ((x.bit_count() & 1) << (L-1)))
        for x in range(G)
    ]
    d = max((j for p in raw for i, j in p), default=0)
    qdeg = max((i for p in raw for i, j in p), default=0)
    # Each item is (Bernstein index j, q-degree k, Fourier row, table row).
    out = []
    for j in range(d+1):
        for k in range(qdeg+1):
            ft = [
                sum(Q(comb(j, r), comb(d, r))*p.get((k, r), 0)
                    for r in range(j+1))
                for p in raw
            ]
            # Full Walsh normalization is 2G; the quotient group has size G.
            table = tuple(Q(x, 2*G) for x in wht(ft))
            out.append((j, k, tuple(ft), table))
    return d, qdeg, out

def dihedral(w):
    n = len(w)
    rev = w[::-1]
    return min(
        [w[i:]+w[:i] for i in range(n)] +
        [rev[i:]+rev[:i] for i in range(n)]
    )

def even_words(L):
    return [w for w in product((1, 2, 3), repeat=L) if sum(w) % 2 == 0]

def screen():
    totals = [0, 0, 0, 0]
    for L in range(2, 7):
        words = even_words(L)
        orbits = {}
        for w in words:
            orbits.setdefault(dihedral(w), []).append(w)
        bern_tables = qvectors = 0
        for rep, members in sorted(orbits.items()):
            d, qdeg, data = bernstein_data(rep)
            bern_tables += len(members)*(d+1)
            qvectors += len(members)*(d+1)*(qdeg+1)
            for j, k, ft, table in data:
                assert min(ft) >= 0
                assert min(table) >= 0
        totals[0] += len(words)
        totals[1] += len(orbits)
        totals[2] += bern_tables
        totals[3] += qvectors
        print("L", L, "ordered", len(words), "dihedral_reps", len(orbits),
              "Bernstein_tables", bern_tables, "q_vectors", qvectors,
              flush=True)
    print("TOTAL", *totals)

def subspaces(d):
    for k in range(d+1):
        for piv in combinations(range(d), k):
            free = [(i, j) for i, p in enumerate(piv)
                    for j in range(p+1, d) if j not in piv]
            for mask in range(1 << len(free)):
                basis = [1 << p for p in piv]
                for h, (i, j) in enumerate(free):
                    if (mask >> h) & 1:
                        basis[i] |= 1 << j
                H = [0]
                for b in basis:
                    H += [x ^ b for x in H]
                yield tuple(H)

def subgroup_columns(d):
    G = 1 << d
    return sorted(set(
        tuple(int(x in H) for x in range(G)) for H in subspaces(d)
    ))

def exact_certificate(cols, target):
    if not any(target):
        return ()
    A = [[float(col[r]) for col in cols] for r in range(len(target))]
    lp = linprog([1.0]*len(cols), A_eq=A,
                 b_eq=[float(x) for x in target],
                 bounds=(0, None), method="highs")
    if not lp.success:
        return None
    active = [i for i, x in enumerate(lp.x) if x > 1e-8]
    M = Matrix([[cols[i][r] for i in active]
                for r in range(len(target))])
    b = Matrix([Rational(x.numerator, x.denominator) for x in target])
    try:
        sol, params = M.gauss_jordan_solve(b)
    except ValueError:
        return None
    if params.rows:
        return None
    coeff = [Q(int(x.p), int(x.q)) for x in sol]
    if any(x < 0 for x in coeff):
        return None
    assert all(
        sum(coeff[t]*cols[i][r] for t, i in enumerate(active)) == target[r]
        for r in range(len(target))
    )
    return [(active[t], coeff[t]) for t in range(len(active)) if coeff[t]]

def certify_reps(lengths, lo6=0, hi6=52):
    cache = {}
    for L in lengths:
        reps = sorted({dihedral(w) for w in even_words(L)})
        if L == 6:
            reps = reps[lo6:hi6]
        cols = subgroup_columns(L-1)
        G = 1 << (L-1)
        vectors = zeros = terms_max = 0
        for labels in reps:
            _, _, data = bernstein_data(labels)
            for j, k, ft, table in data:
                vectors += 1
                if not any(table):
                    zeros += 1
                    continue
                key = (L, table)
                if key not in cache:
                    cache[key] = exact_certificate(cols, table)
                cert = cache[key]
                assert cert is not None, (labels, j, k)
                terms_max = max(terms_max, len(cert))
        print("CERT", L, "reps", len(reps), "vectors", vectors,
              "zero", zeros, "max_terms", terms_max, flush=True)

def subspaces_inside(support):
    out = {frozenset({0})}
    for v in sorted(support - {0}):
        for H in tuple(out):
            H2 = H | frozenset(x ^ v for x in H)
            if H2 <= support:
                out.add(H2)
    return sorted(out, key=lambda H: (len(H), tuple(sorted(H))))

def exact_cert_supported(g, Hs, extra=None):
    G = len(g)
    cols = [tuple(int(x in H) for x in range(G)) for H in Hs]
    if extra is not None:
        cols.append(tuple(extra))
    supp = {x for x, v in enumerate(g) if v}
    cols = [c for c in cols if all(c[x] == 0 for x in range(G)
                                   if x not in supp)]
    return exact_certificate(cols, g)

def boundary(labels, add_sphere=False):
    d, qdeg, data = bernstein_data(labels)
    L = len(labels)
    G = 1 << (L-1)
    failures = repairs = 0
    terms_max = 0
    p1 = [int(x.bit_count() == 1) for x in range(G)]
    p1conv = tuple(sum(p1[y]*p1[x ^ y] for y in range(G))
                   for x in range(G))
    for j, k, ft, table in data:
        if not any(table):
            continue
        support = {x for x, z in enumerate(table) if z}
        Hs = subspaces_inside(support)
        cert = exact_cert_supported(table, Hs)
        if cert is None:
            failures += 1
            sphere = p1conv if all(p1conv[x] == 0
                                   for x in range(G) if x not in support) else None
            assert add_sphere and sphere is not None
            cert = exact_cert_supported(table, Hs, sphere)
            assert cert is not None
            repairs += 1
        terms_max = max(terms_max, len(cert))
    print("BOUNDARY", labels, "d", d, "qdeg", qdeg,
          "vectors", len(data), "subgroup_misses", failures,
          "sphere_repairs", repairs, "max_terms", terms_max, flush=True)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "screen"
    if mode in ("-h", "--help"):
        print("modes: screen | smallcert | sixcert LO HI | boundary7 | boundary9")
    elif mode == "screen":
        screen()
    elif mode == "smallcert":
        certify_reps(range(2, 6))
    elif mode == "sixcert":
        certify_reps((6,), int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "boundary7":
        boundary((1,1,1,1,2,2,2))
    elif mode == "boundary9":
        boundary((1,)*8 + (6,), add_sphere=True)
    else:
        raise SystemExit("unknown mode")
