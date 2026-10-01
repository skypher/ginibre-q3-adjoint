from collections import defaultdict

def consumer_value(word):
    # State (left spin, right spin) with exact integer multiplicities.
    states = {(0, 0): 1}
    for signed_n in word:
        n = abs(signed_n)
        sign = 1 if signed_n > 0 else -1
        nxt = defaultdict(int)
        for (x, y), value in states.items():
            for spin in range(abs(x - n), x + n + 1, 2):
                nxt[spin, y] += value
            for spin in range(abs(y - n), y + n + 1, 2):
                nxt[x, spin] += sign * value
        states = nxt
    return states.get((0, 0), 0)

def distance(word):
    labels = [abs(x) for x in word]
    return (sum(labels) - 2 * max(labels)) // 2

smallest = [
    ([-1] * 14 + [1, 2, -4, -5], 3952568),
    ([-1] * 14 + [1, -3, -3, 5], 6028200),
    ([-1] * 14 + [1, -3, 3, -5], 2789224),
    ([-1] * 14 + [-3, -4, 5], 16211880),
    ([-1] * 14 + [3, -4, -5], 6731864),
]
for word, value in smallest:
    assert sum(x < 0 for x in word) % 2 == 0
    assert sum(map(abs, word)) == 26 and distance(word) == 8
    assert consumer_value(word) == value

shapes = {
    8: [
        (2242, [-1] * 22 + [-8, -8, 22], 47767619552),
        (2234, [-1] * 22 + [1, -7, -8, 22], 15651154908),
        (2149, [-1] * 24 + [-6, -8, 22], 225638266660),
    ],
    9: [
        (3124, [-1] * 22 + [-8, -9, 21], 129046630196),
        (3111, [-1] * 22 + [1, -7, -9, 21], 34296133520),
        (3018, [-1] * 20 + [1, -9, -9, 21], 4407933090),
    ],
    10: [
        (4175, [-1] * 22 + [-8, -10, 20], 263307351560),
        (4164, [-1] * 22 + [1, -7, -10, 20], 58243687612),
        (4066, [-1] * 24 + [-6, -10, 20], 1894253713990),
    ],
}
for delta, entries in shapes.items():
    for frequency, word, value in entries:
        assert sum(x < 0 for x in word) % 2 == 0
        assert sum(map(abs, word)) == 60 and distance(word) == delta
        assert consumer_value(word) == value

def check_hist(start, counts, total):
    assert sum(counts) == total
    return dict(zip(range(start, start + len(counts)), counts))

phase0_d8 = {
    "h": {0: 29190},
    "D": check_hist(21, [4990, 9372, 14828], 29190),
    "L": check_hist(6, [5, 389, 1518, 2819, 3915, 4357, 4252, 3730,
                        3035, 2234, 1464, 842, 409, 163, 50, 6, 2], 29190),
    "min_large": {0: 29190},
}
phase1_d8 = {
    "h": {0: 202835, 1: 72325, 2: 48},
    "D": check_hist(16, [
        2140, 3022, 3412, 4172, 4674, 6830, 8419, 9984, 11631, 12641,
        13560, 13716, 13588, 12372, 10790, 11959, 13685, 15029, 16986,
        18543, 20687, 22460, 24908], 275208),
    "L": check_hist(8, [
        33, 328, 1106, 2799, 4738, 7339, 10072, 13153, 15738, 18407,
        20292, 21844, 22292, 22184, 20923, 18806, 16393, 14120, 11763,
        9639, 7560, 5776, 4135, 2779, 1646, 842, 352, 117, 22, 10], 275208),
    "min_large": dict(zip(
        [0] + list(range(9, 23)),
        [11740, 23695, 20787, 17614, 16240, 15109, 15073, 15107,
         16031, 16808, 18316, 19495, 21310, 22826, 25057])),
}
assert all(sum(hist.values()) == 275208 for hist in phase1_d8.values())

structures = {
    0: {8: (1002, 28188, 29190),
        9: (881, 27855, 28736),
        10: (517, 17116, 17633)},
    1: {8: (16052, 259156, 275208),
        9: (22480, 336145, 358625),
        10: (29534, 422193, 451727)},
}
for rows in structures.values():
    for two_core, mech69_remainder, total in rows.values():
        assert two_core + mech69_remainder == total

phase0_marginal = [23706, 292, 0, 1186, 0]
phase1_marginal = [204738, 27147, 0, 8041, 0]
assert 100743 - sum(phase0_marginal) == 75559
assert 9823286 - sum(phase1_marginal) == 9583360
assert sum(phase0_marginal) == 25184
assert sum(phase1_marginal) == 239926

print("EXACT AUDIT PASS")