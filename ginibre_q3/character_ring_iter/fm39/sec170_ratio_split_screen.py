from collections import defaultdict
from fractions import Fraction

def coefficient(word, a, b):
    # [U_a(x) U_b(y)] product_i (U_|n_i|(x) + sign(n_i) U_|n_i|(y)).
    state = {(0, 0): 1}
    for z in word:
        n, eps = abs(z), (1 if z > 0 else -1)
        nxt = defaultdict(int)
        for (i, j), v in state.items():
            for k in range(abs(i-n), i+n+1, 2):
                nxt[(k, j)] += v
            for k in range(abs(j-n), j+n+1, 2):
                nxt[(i, k)] += eps*v
        state = nxt
    return state.get((a, b), 0)

def g(B, p):
    return coefficient(B, abs(p), 0)

def top_pair(B):
    candidates = []
    for i in range(len(B)):
        for j in range(i+1, len(B)):
            a, b = B[i], B[j]
            if (abs(a)+abs(b)) % 2 == 0:
                candidates.append(((abs(a)+abs(b), max(abs(a), abs(b))), a, b))
    _, a, b = max(candidates)
    return a, b

def remove_pair(B, pair):
    C = list(B)
    C.remove(pair[0])
    C.remove(pair[1])
    return tuple(C)

def check_residual(B, p):
    W = sum(map(abs, B))
    mx = max(map(abs, B))
    delta = (W-abs(p)) // 2
    signs = {}
    for z in B:
        signs.setdefault(abs(z), set()).add(1 if z > 0 else -1)
    assert all(len(v) == 1 for v in signs.values())
    sigma = -1 if sum(z < 0 for z in B) % 2 else 1
    assert p == sigma*abs(p)
    assert W % 2 == abs(p) % 2 and W > abs(p)
    assert abs(p) >= max(6, mx) and delta >= 8 and mx <= delta
    assert sum(abs(z) >= 3 for z in B) >= 2
    return W, delta

def flip_values(B, p):
    factors = list(B) + [p]
    values = []
    for z in factors:
        if z not in values:
            values.append(z)
    out = []
    for i, a in enumerate(values):
        for b in values[i:]:
            if a == b and factors.count(a) < 2:
                continue
            C = factors.copy()
            C.remove(a)
            C.remove(b)
            D = coefficient(C, abs(a), abs(b)) * (1 if b > 0 else -1)
            out.append((a, b, D))
    return out

B43 = tuple(map(int,
    "1 2 2 2 1 1 2 1 2 2 2 1 1 -3 2 1 1 1 2 1 1 2 -3 2 2 1 1".split()))
p43 = 7
W, delta = check_residual(B43, p43)
pair = top_pair(B43)
assert (W, delta, pair) == (43, 18, (-3, -3))
assert Fraction(sum(map(abs, pair)), delta) == Fraction(1, 3)
assert (g(B43, p43), g(remove_pair(B43, pair), p43)) == (
    28362093763816, 29753344469664)
assert (1, 1, 14086016005744) in flip_values(B43, p43)
print("W43:", W, delta, pair, Fraction(1, 3),
      g(B43, p43), g(remove_pair(B43, pair), p43),
      "D(1,1)=14086016005744")

B111 = tuple(map(int,
    "3 2 3 2 2 1 2 4 3 3 3 4 2 3 -5 3 -5 3 3 4 4 3 4 3 4 4 3 4 2 3 4 4 2 3 4".split()))
p111 = 9
W, delta = check_residual(B111, p111)
pair = top_pair(B111)
assert (W, delta, pair) == (111, 51, (-5, -5))
assert Fraction(10, 51) < Fraction(1, 2)
assert g(B111, p111) == 4733800780161150966287804
assert g(remove_pair(B111, pair), p111) == 4827798007122625669517662
assert (2, 1, 2364673630242973617226396) in flip_values(B111, p111)
print("W111:", W, delta, pair, Fraction(10, 51),
      g(B111, p111), g(remove_pair(B111, pair), p111),
      "D(2,1)=2364673630242973617226396")

B48 = tuple(map(int,
    "-1 -1 -2 -3 -3 -3 -3 -4 -4 -4 -5 -5 -5 -5".split()))
p48 = 8
W, delta = check_residual(B48, p48)
pair = top_pair(B48)
D = flip_values(B48, p48)
assert (W, delta, pair) == (48, 20, (-5, -5))
assert Fraction(sum(map(abs, pair)), delta) == Fraction(1, 2)
assert g(B48, p48) == 4075371
assert g(remove_pair(B48, pair), p48) == 160741
assert len(D) == 19 and max(v for _, _, v in D) == -20406
print("W48:", W, delta, pair, Fraction(1, 2),
      g(B48, p48), g(remove_pair(B48, pair), p48),
      "flip_pairs=19", "max_D=-20406")
