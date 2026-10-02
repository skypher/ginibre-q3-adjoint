from collections import defaultdict
from fractions import Fraction

def coefficient(word, a, b):
    """Coefficient of U_a(x) U_b(y) in the signed character product."""
    word = sorted(word, key=lambda z: -abs(z))
    remaining = sum(abs(z) for z in word)
    state = {(0, 0): 1}

    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        remaining -= n
        nxt = defaultdict(int)

        for (s, t), value in state.items():
            if abs(t - b) <= remaining:
                for s2 in range(abs(s - n), s + n + 1, 2):
                    if abs(s2 - a) <= remaining:
                        nxt[(s2, t)] += value

            if abs(s - a) <= remaining:
                for t2 in range(abs(t - n), t + n + 1, 2):
                    if abs(t2 - b) <= remaining:
                        nxt[(s, t2)] += eps * value

        state = nxt

    return state.get((a, b), 0)

def g(B, p):
    return coefficient(B, abs(p), 0)

def top_pair(B):
    choices = [
        (abs(a) + abs(b), max(abs(a), abs(b)), a, b)
        for i, a in enumerate(B)
        for b in B[i + 1:]
        if (abs(a) + abs(b)) % 2 == 0
    ]
    _, _, a, b = max(choices)
    return a, b

B = (-40, 42, -44, 46, 48, 50, 52, 54, 56, 58, 60)
p = 62
W = sum(map(abs, B))
delta = (W - abs(p)) // 2
pair = top_pair(B)

assert W == 550
assert delta == 244
assert pair == (58, 60)
assert Fraction(sum(map(abs, pair)), delta) == Fraction(59, 122)

child = list(B)
child.remove(pair[0])
child.remove(pair[1])
parent_value = g(B, p)
child_value = g(child, p)

assert parent_value == 453207534222864
assert child_value == 179646349605

factors = list(B) + [p]
differences = []
for i, a in enumerate(factors):
    for j in range(i + 1, len(factors)):
        b = factors[j]
        rest = factors[:i] + factors[i + 1:j] + factors[j + 1:]
        d = coefficient(rest, abs(a), abs(b))
        if b < 0:
            d = -d
        differences.append((a, b, d))
        print("pair", len(differences), "/ 66", a, b,
              "D =", d, flush=True)

assert len(differences) == 66
max_D = max(d for _, _, d in differences)
assert max_D == -34349665
assert all(d < 0 for _, _, d in differences)

print("W =", W, "delta =", delta, "TopPair =", pair)
print("wTP/delta =", Fraction(sum(map(abs, pair)), delta))
print("g_p(B) =", parent_value)
print("g_p(B-TP) =", child_value)
print("g_p(B-TP)/g_p(B) =", Fraction(child_value, parent_value))
print("flip pairs =", len(differences), "max D =", max_D)
