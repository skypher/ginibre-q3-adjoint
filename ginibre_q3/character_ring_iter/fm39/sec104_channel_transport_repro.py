"""FM-SEC104 (luna_max_vesta): equal-label singlet channel transports E_- -> E_+; rank losses at (1,1,1,1)."""
from functools import lru_cache
from itertools import product
import sympy as sp

@lru_cache(None)
def paths(labels):
    """Fusion paths for SU(2), ending at highest weight 0."""
    states = {0: [()]}
    for n in labels:
        nxt = {}
        for a, old_paths in states.items():
            for c in range(abs(a - n), a + n + 1, 2):
                nxt.setdefault(c, []).extend(p + (c,) for p in old_paths)
        states = nxt
    return tuple(states.get(0, ()))

def side_paths(ns, mask):
    L = len(ns)
    A = tuple(ns[i] for i in range(L) if mask >> i & 1)
    B = tuple(ns[i] for i in range(L) if not (mask >> i & 1))
    return paths(A), paths(B)

def canonical(mask, L):
    return mask if not (mask >> (L - 1) & 1) else mask ^ ((1 << L) - 1)

def signed_slots(ns, T):
    L = len(ns)
    reps = range(1 << (L - 1))
    dims = {s: len(side_paths(ns, s)[0]) * len(side_paths(ns, s)[1])
            for s in reps}
    neg = [s for s in reps if (s & T).bit_count() % 2 and dims[s]]
    pos = [s for s in reps if not ((s & T).bit_count() % 2) and dims[s]]
    return pos, neg, dims

def transport_matrix(ns, T, rule):
    """Exact matrix in external fusion-path coordinates."""
    L = len(ns)
    pos, neg, dims = signed_slots(ns, T)
    row0, col0 = {}, {}
    q = 0
    for s in pos:
        row0[s] = q
        q += dims[s]
    nrow = q
    q = 0
    for s in neg:
        col0[s] = q
        q += dims[s]
    M = sp.zeros(nrow, q)

    for s in neg:
        edges = []
        for i in range(L):
            if not (T >> i & 1):
                continue
            for j in range(L):
                if (T >> j & 1) or ns[i] != ns[j]:
                    continue
                if bool(s >> i & 1) == bool(s >> j & 1):
                    continue
                t = canonical(s ^ (1 << i) ^ (1 << j), L)
                if t in row0:
                    edges.append((i, j, t))
        if rule == "first" and edges:
            edges = [min(edges, key=lambda e: (e[0], e[1]))]

        srcA, srcB = side_paths(ns, s)
        src_basis = [(a, b) for a in srcA for b in srcB]
        for i, j, t in edges:
            dstA, dstB = side_paths(ns, t)
            dst_basis = [(a, b) for a in dstA for b in dstB]
            raw = s ^ (1 << i) ^ (1 << j)
            complement_flip = raw != t
            coeff = (-1) ** ns[i] if rule == "FS" else 1
            for k, (a, b) in enumerate(src_basis):
                target = (b, a) if complement_flip else (a, b)
                M[row0[t] + dst_basis.index(target), col0[s] + k] += coeff
    return pos, neg, M

for ns, T in [((1,1,1,1), 3),
              ((1,1,1,1,2), 3),
              ((1,5,2,2), 3)]:
    pos, neg, dims = signed_slots(ns, T)
    dplus = sum(dims[s] for s in pos)
    dminus = sum(dims[s] for s in neg)
    print("boundary", ns, "T=", T, "E+/E-=", dplus, dminus)
    for rule in ("sum", "FS", "first"):
        _, _, M = transport_matrix(ns, T, rule)
        print(rule, "rank=", M.rank(), "kernel=", M.nullspace())

for L in range(1, 4):
    for ns in product((1, 2, 3), repeat=L):
        if sum(ns) % 2:
            continue
        for T in range(1 << L):
            if T.bit_count() % 2:
                continue
            for s in range(1 << (L - 1)):
                if (s & T).bit_count() % 2:
                    A = tuple(ns[i] for i in range(L) if s >> i & 1)
                    B = tuple(ns[i] for i in range(L) if not (s >> i & 1))
                    assert len(paths(A)) * len(paths(B)) == 0
print("verified E-=0 for every ordered list with L<=3, labels<=3")
