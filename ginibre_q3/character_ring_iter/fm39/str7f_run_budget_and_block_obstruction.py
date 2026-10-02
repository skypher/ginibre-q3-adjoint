from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
import sys

if any(a in ("-h", "--help") for a in sys.argv[1:]):
    print("Exact FM-STR7f verifier: checks the bounded-d Sp(4) cut obstruction and small residual-run budgets.")
    raise SystemExit(0)

def stamp(msg):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), msg, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def add(dst, src, scale=1):
    for key, val in src.items():
        dst[key] += scale*val

def mul_basis(f, a, b):
    out = defaultdict(int)
    for (u, v), c in f.items():
        for x in cg(u, a):
            for y in cg(v, b):
                out[x, y] += c
    return {k:v for k,v in out.items() if v}

def mul_d(f):
    out = defaultdict(int)
    add(out, mul_basis(f, 1, 0))
    add(out, mul_basis(f, 0, 1), -1)
    return {k:v for k,v in out.items() if v}

def mul_K(f, n):
    out = defaultdict(int)
    for j in range(n):
        add(out, mul_basis(f, j, n-1-j))
    return {k:v for k,v in out.items() if v}

def mul_label(f, signed_n):
    n = abs(signed_n)
    out = defaultdict(int)
    add(out, mul_basis(f, n, 0))
    add(out, mul_basis(f, 0, n), 1 if signed_n > 0 else -1)
    return {k:v for k,v in out.items() if v}

def row(word):
    f = {(0,0):1}
    for z in sorted(word, key=abs, reverse=True):
        f = mul_label(f, z)
    return f

def sp4_char(A, B):
    out = defaultdict(int)
    for j in range(B+1):
        for k in range(A-B+1):
            out[j+k, j+A-B-k] += 1
    return out

def sp4_expand(f):
    f = defaultdict(int, f)
    assert all(f.get((u,v),0) == f.get((v,u),0)
               for u,v in f), "not swap-symmetric"
    out = defaultdict(int)
    while any(f.values()):
        D = max(u+v for (u,v),c in f.items() if c)
        top = [f.get((D-b,b),0) for b in range(D//2+1)]
        prev = 0
        for b, val in enumerate(top):
            mult = val-prev
            prev = val
            if mult:
                A = D-b
                out[A,b] += mult
                for key, coeff in sp4_char(A,b).items():
                    f[key] -= mult*coeff
        f = defaultdict(int, {k:v for k,v in f.items() if v})
    return {k:v for k,v in out.items() if v}

def block_expansion(labels, residual_d_power):
    f = {(0,0):1}
    for _ in range(residual_d_power):
        f = mul_d(f)
    for n in labels:
        f = mul_K(f, n)
    return sp4_expand(f)

def bounded_d_lift_exists(labels):
    r = len(labels)
    for extracted in range(3):
        q = r-extracted
        if q < 0 or q % 2:
            continue
        if all(c >= 0 for c in block_expansion(labels, q).values()):
            return True
    return False

stamp("begin exact Sp(4) cut check; Python integers")
full = tuple(-n for n in range(1,9))
assert row(full).get((0,0),0) == 980
B = full[:-1]
p = 8
assert len(B) == 7 and sum(abs(z) for z in B) == 28
assert (-1)**sum(z < 0 for z in B) * p == -8
assert (28-p) % 2 == 0 and (28-p)//2 == 10
assert p >= max(abs(z) for z in B) and max(abs(z) for z in B) <= 10
assert sum(abs(z) >= 3 for z in B) >= 2
Ds = []
for i,j in combinations(range(8),2):
    C = [full[t] for t in range(8) if t not in (i,j)]
    Ds.append(-row(C).get((i+1,j+1),0))
assert len(Ds) == 28 and min(Ds) == -65 and max(Ds) == -8
assert all(D < 0 for D in Ds)
stamp("residual (-1,...,-8): Phi=980; all 28 flips lie in [-65,-8]")

checked = 0
for r in range(3,7):
    for labels in combinations(range(1,9), r):
        admissible = [e for e in range(3)
                      if r >= e and (r-e)%2 == 0]
        for e in admissible:
            coeffs = block_expansion(labels, r-e)
            assert any(c < 0 for c in coeffs.values()), (labels,e,coeffs)
            checked += 1
        assert not bounded_d_lift_exists(labels)
    stamp(f"checked all size-{r} all-minus blocks from labels 1..8")
assert checked == 308
assert block_expansion((1,2,3),2) == {
    (5,0):1, (4,1):-1, (3,2):-1, (3,0):2, (1,0):2
}
assert block_expansion((1,3,8),2).get((10,1)) == -1
assert block_expansion((2,4,6),2).get((9,2)) == -1
stamp("210 blocks; 308 admissible d^e quotients virtual; TopPair blocks have chi_(10,1), chi_(9,2) coefficients -1")

assignments = 0
for i,j in combinations(range(8),2):
    remaining = tuple(n for t,n in enumerate(range(1,9)) if t not in (i,j))
    for mask in range(1 << 6):
        left = tuple(remaining[q] for q in range(6) if mask & (1 << q))
        right = tuple(remaining[q] for q in range(6) if not (mask & (1 << q)))
        large = left if len(left) >= len(right) else right
        assert len(large) >= 3 and not bounded_d_lift_exists(large)
        assignments += 1
assert assignments == 1792
stamp("all 28 pair deletions and 1,792 two-block assignments fail this cut class")

def split(C):
    A, B = [], []
    wa = wb = 0
    for z in sorted(C, key=abs, reverse=True):
        if wa <= wb:
            A.append(z); wa += abs(z)
        else:
            B.append(z); wb += abs(z)
    if wa > wb:
        A, B = B, A
    return A, B

def run_budget(word, i, j):
    C = [z for k,z in enumerate(word) if k not in (i,j)]
    A, B = split(C)
    X = row(A)
    Y = row(B + [word[i],word[j]])
    Z = row(B + [-word[i],-word[j]])
    max_t = max((u+v for u,v in X), default=0)
    P = [0]*(max_t+1)
    M = [0]*(max_t+1)
    for (u,v), x in X.items():
        q, f = Y.get((u,v),0), Z.get((u,v),0)
        assert (q+f)%2 == 0 and (q-f)%2 == 0
        P[u+v] += x*((q+f)//2)
        M[u+v] += x*((q-f)//2)
    prefix = 0
    vals = []
    for a,b in zip(P,M):
        prefix += a + min(b,0)
        vals.append(prefix)
    return min(vals, default=0), sum(P)+sum(M), A, B, P, M, vals

run_cases = 0
layer_checks = 0
run_receipt = None
for k in range(7,11):
    W = k*(k+1)//2
    sigma = -1 if k%2 else 1
    ps = [q for q in range(max(k+1,6),k+6)
          if (W-q)%2 == 0 and (W-q)//2 >= 8 and k <= (W-q)//2]
    for q in ps:
        word = tuple(-n for n in range(1,k+1)) + (sigma*q,)
        i,j = k-3,k-1
        assert (abs(word[i])+abs(word[j]))%2 == 0
        bt, phi, A, C, P, M, vals = run_budget(word,i,j)
        assert bt >= 0 and phi == row(word).get((0,0),0)
        for pt,mt in zip(P,M):
            assert pt >= 0 and pt+mt >= 0
            layer_checks += 1
        run_cases += 1
        if k == 7 and q == 8:
            assert word == full and (A,C) == ([-8,-3,-1],[-6,-4,-2])
            assert phi == 980
            assert vals == [0,0,0,0,118,118,394,394,674,674,888,888,980]
            run_receipt = vals
    stamp(f"residual run checks through k={k}; admissible p values={ps}")
assert run_cases == 10 and run_receipt is not None
stamp(f"run layerwise checks={layer_checks}; k=7,p=8 B_T={run_receipt}")
stamp(f"FM-STR7f verifier PASS; exact residual run cases={run_cases}; no uniform proof claimed")
