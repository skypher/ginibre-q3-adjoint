# FM-SEC68 (luna_max_eris) printed code: Hypothesis B fails at (1^8, 6).
from functools import lru_cache

@lru_cache(None)
def fusion(labels):
    b = {0: 1}
    for n in labels:
        out = {}
        for p, v in b.items():
            for q in range(abs(p - n), p + n + 1, 2):
                out[q] = out.get(q, 0) + v
        b = out
    return b

def inv(labels):
    return fusion(tuple(labels)).get(0, 0)

def table_one_high(a):
    labels = (1,) * 8 + (a,)
    f = []
    for x in range(1 << 8):
        left = tuple(labels[i] for i in range(8) if x >> i & 1)
        right = tuple(labels[i] for i in range(8)
                      if not (x >> i & 1)) + (a,)
        f.append(inv(left) * inv(right))
    return f

def subspaces(d):
    layer = {(0,)}
    yield (0,)
    for _ in range(d):
        nxt = set()
        for H in layer:
            seen = set(H)
            for v in range(1, 1 << d):
                if v in seen:
                    continue
                coset = {v ^ x for x in H}
                seen.update(coset)
                nxt.add(tuple(sorted(set(H) | coset)))
        for H in sorted(nxt):
            yield H
        layer = nxt

for a in (5, 7):
    assert not any(table_one_high(a))
    print("length 9, (1^8,%d): zero table, empty certificate" % a)

f8 = table_one_high(8)
assert f8 == [1] + [0] * 255
print("length 9, (1^8,8): certificate 1_{0}, coefficient 1")

f = table_one_high(6)
assert all(f[x] == (7 if x == 0 else 1 if x.bit_count() == 2 else 0)
           for x in range(1 << 8))

y_by_weight = (3, -1, -1, 3, 3, -1, -1, 3, 3)
y = lambda x: y_by_weight[x.bit_count()]
separation = sum(y(x) * f[x] for x in range(1 << 8))
assert separation == -7

F = f[:]
h = 1
while h < len(F):
    for i in range(0, len(F), 2 * h):
        for j in range(h):
            a0, b0 = F[i + j], F[i + h + j]
            F[i + j], F[i + h + j] = a0 + b0, a0 - b0
    h *= 2
assert min(F) == 3

count = zero_slacks = 0
minimum = None
for count, H in enumerate(subspaces(8), 1):
    slack = sum(y(x) for x in H)
    minimum = slack if minimum is None else min(minimum, slack)
    zero_slacks += (slack == 0)
assert count == 417199 and minimum == 0 and zero_slacks == 11830

print("subspaces checked:", count)
print("minimum subspace sum:", minimum)
print("sum_x y(x) f(x):", separation)
print("Fourier minimum:", min(F))
print("zero-sum subspaces:", zero_slacks)
