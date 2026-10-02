from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction

def progress(message):
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    print(stamp, message, flush=True)

def fuse(a, b):
    return range(abs(a - b), a + b + 1, 2)

def multiply_factor(table, n, sign):
    out = defaultdict(int)
    for (r, s), coefficient in table.items():
        for t in fuse(r, n):
            out[t, s] += coefficient
        for t in fuse(s, n):
            out[r, t] += sign * coefficient
    return {key: value for key, value in out.items() if value}

def coefficients(word):
    table = {(0, 0): 1}
    for n, sign in word:
        table = multiply_factor(table, n, sign)
    return table

def pair_products(left, right):
    a = coefficients(left)
    b = coefficients(right)
    return {key: value * b.get(key, 0) for key, value in a.items()
            if value * b.get(key, 0)}

def height_profile(channels):
    out = defaultdict(int)
    for (r, s), value in channels.items():
        out[r + s] += value
    return dict(out)

def heat_profile(channels):
    out = defaultdict(int)
    for (r, s), value in channels.items():
        out[r * (r + 2) + s * (s + 2)] += value
    return dict(out)

def evaluate(profile, q):
    return sum(Fraction(value) * q**height
               for height, value in profile.items())

def derivative_q(profile, q):
    return sum(Fraction(height * value) * q**(height - 1)
               for height, value in profile.items() if height)

def abel_reconstruction(profile, q):
    top = max(profile, default=0)
    cumulative = 0
    result = Fraction(0)
    for height in range(top):
        cumulative += profile.get(height, 0)
        result += (1 - q) * q**height * cumulative
    return result + q**top * sum(profile.values())

def pair_free(word):
    signs = {}
    for n, sign in word:
        if n in signs and signs[n] != sign:
            return False
        signs[n] = sign
    return sum(sign < 0 for _, sign in word) % 2 == 0

progress('start exact SU(2) fusion checks')
minus_one_block = coefficients(((1, -1), (1, -1)))
assert minus_one_block == {
    (0, 0): 2, (2, 0): 1, (0, 2): 1, (1, 1): -2
}

# Pair-free list (-3,-3), split into singleton halves.
a1 = ((3, -1),)
b1 = ((3, -1),)
ch1 = pair_products(a1, b1)
p1 = height_profile(ch1)
assert pair_free(a1 + b1)
assert ch1 == {(3, 0): 1, (0, 3): 1}
assert p1 == {3: 2}
assert evaluate(p1, Fraction(1)) == 2
assert derivative_q(p1, Fraction(1)) == 6
assert coefficients(a1).get((0, 0), 0) * coefficients(b1).get((0, 0), 0) == 0

# Pair-free all-minus list (-3)^4 (-1)^2 (-2)^2, split by label classes.
a2 = ((3, -1),) * 4
b2 = ((1, -1),) * 2 + ((2, -1),) * 2
ch2 = pair_products(a2, b2)
p2 = height_profile(ch2)
assert pair_free(a2 + b2)
assert ch2[(3, 3)] == -128
assert p2 == {0: 84, 2: 120, 4: 100, 6: -108}
assert [sum(p2.get(h, 0) for h in range(T + 1))
        for T in range(7)] == [84, 84, 204, 204, 304, 304, 196]
assert coefficients(a2).get((0, 0), 0) * coefficients(b2).get((0, 0), 0) == 84
assert evaluate(p2, Fraction(1)) == 196
assert derivative_q(p2, Fraction(1)) == -8
for q in (Fraction(0), Fraction(1, 4), Fraction(1, 2),
          Fraction(3, 4), Fraction(1)):
    assert evaluate(p2, q) == abel_reconstruction(p2, q)
    assert evaluate(p2, q) >= 0
progress('virtual block, pair-free Poisson profiles, and Abel identity checked')

# Pair-free heat examples with opposite exact t-derivative signs at t=0.
h1 = heat_profile(ch1)
assert h1 == {15: 2}
dt0_1 = -sum(E * value for E, value in h1.items())
assert dt0_1 == -30

a3 = ((3, 1),) * 2
b3 = ((1, -1),) * 2 + ((2, 1),) * 4
ch3 = pair_products(a3, b3)
h3 = heat_profile(ch3)
assert pair_free(a3 + b3)
assert h3 == {0: 112, 8: 204, 24: 138, 48: 52, 30: -256}
dt0_3 = -sum(E * value for E, value in h3.items())
assert dt0_3 == 240
progress('opposite heat-derivative signs checked on pair-free lists')
print('q profiles:', p1, p2, flush=True)
print('q derivatives at 1:', derivative_q(p1, Fraction(1)),
      derivative_q(p2, Fraction(1)), flush=True)
print('heat profiles by Casimir:', h1, h3, flush=True)
print('heat derivatives at t=0:', dt0_1, dt0_3, flush=True)
progress('PASS exact verifier')
