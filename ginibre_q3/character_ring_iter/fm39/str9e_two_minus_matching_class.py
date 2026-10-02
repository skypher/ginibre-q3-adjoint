import sys
if any(arg in ("-h", "--help") for arg in sys.argv[1:]):
    print("FM-STR9e verifier: Sp(4) decomposition of K_n S_a and two-minus interior-prefix checks.")
    raise SystemExit(0)
from datetime import datetime, timezone
from itertools import combinations_with_replacement

def log(*xs):
    print(datetime.now(timezone.utc).isoformat(timespec="seconds"), *xs, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def multiply_factor(table, n, eps=1):
    out = {}
    for (r, s), value in table.items():
        for q in cg(r, n):
            out[q, s] = out.get((q, s), 0) + value
        for q in cg(s, n):
            out[r, q] = out.get((r, q), 0) + eps * value
    return {key: value for key, value in out.items() if value}

def character(factors):
    table = {(0, 0): 1}
    for label, sign in factors:
        table = multiply_factor(table, label, sign)
    return table

def K_table(n):
    return {(i, n-1-i): 1 for i in range(n)}

def d_times(table):
    return multiply_factor(table, 1, -1)

def expected_KS_decomposition(n, a):
    out = {}
    for j in range(min(n, a) + 1):
        q = n + a - 1 - 2*j
        if q >= 0:
            out[q, 0] = out.get((q, 0), 0) + 1
    if n > a:
        out[n-1, a] = out.get((n-1, a), 0) + 1
    elif a > n:
        out[a-1, n] = out.get((a-1, n), 0) - 1
    return {key: value for key, value in out.items() if value}

def actual_decomposition(table):
    numerator = d_times(table)
    for (i, j), value in numerator.items():
        assert numerator.get((j, i), 0) == -value
    return {(i-1, j): value for (i, j), value in numerator.items()
            if i > j and value}

def check_KS_identity(max_n=18, max_a=20):
    checks = 0
    for n in range(1, max_n+1):
        for a in range(max_a+1):
            lhs = actual_decomposition(multiply_factor(K_table(n), a, 1))
            rhs = expected_KS_decomposition(n, a)
            assert lhs == rhs, (n, a, lhs, rhs)
            checks += 1
        if n % 3 == 0:
            log("K_n S_a decomposition", "n", n, "through a", max_a,
                "checks", checks)
    assert checks == max_n * (max_a+1)
    return checks

def prefix_values(A, Bprime):
    fa, fb = character(A), character(Bprime)
    hmax = max((r+s for r, s in fa), default=0)
    levels = [0] * (hmax+1)
    for (r, s), value in fa.items():
        levels[r+s] += value * fb.get((r, s), 0)
    prefixes = []
    total = 0
    for value in levels:
        total += value
        prefixes.append(total)
    full = character(A+Bprime).get((0, 0), 0)
    assert total == full, (A, Bprime, total, full)
    assert all(value >= 0 for value in prefixes), (A, Bprime, prefixes)
    return prefixes, full

def matched_layout(n, m, exceptional):
    if not exceptional:
        return None
    if len(exceptional) == 1:
        a = exceptional[0]
        if a <= n:
            return n, a, m, None
        if a <= m:
            return m, a, n, None
        return False
    assert len(exceptional) == 2
    a, b = exceptional
    for minus_B, exception_B, minus_A, exception_A in (
        (n, a, m, b), (n, b, m, a),
        (m, a, n, b), (m, b, n, a)):
        if exception_B <= minus_B and exception_A <= minus_A:
            return minus_B, exception_B, minus_A, exception_A
    return False

def verify_prefix_class(max_minus=8, max_exception=8, max_ones=4):
    cases = 0
    min_prefix = None
    for n in range(1, max_minus+1):
        for m in range(1, max_minus+1):
            exceptional_lists = [()]
            exceptional_lists += [(a,) for a in range(2, max_exception+1)]
            exceptional_lists += list(combinations_with_replacement(
                range(2, max_exception+1), 2))
            for exceptional in exceptional_lists:
                layout = matched_layout(n, m, exceptional)
                if layout is False:
                    continue
                for ones in range(max_ones+1):
                    left_ones = ones // 2
                    right_ones = ones-left_ones
                    if layout is None:
                        A = [(1, 1)]*left_ones
                        Bprime = [(n, -1), (m, -1)] + [(1, 1)]*right_ones
                    else:
                        minus_B, exception_B, minus_A, exception_A = layout
                        A = [(minus_A, -1)]
                        if exception_A is not None:
                            A.append((exception_A, 1))
                        A += [(1, 1)]*left_ones
                        Bprime = [(minus_B, -1), (exception_B, 1)]
                        Bprime += [(1, 1)]*right_ones
                    prefixes, _ = prefix_values(A, Bprime)
                    local_min = min(prefixes)
                    min_prefix = local_min if min_prefix is None else min(min_prefix, local_min)
                    cases += 1
        log("prefix class", "minus n", n, "cumulative exact cases", cases)
    return cases, min_prefix

def main():
    log("start exact FM-STR9e verifier")
    identities = check_KS_identity()
    log("PASS K_n S_a identity", identities, "exact decomposition checks")

    # At n=2, the only multiset of two labels >1 with each label <= n is (2,2).
    smallest_bounded = K_table(2)
    smallest_bounded = multiply_factor(smallest_bounded, 2, 1)
    smallest_bounded = multiply_factor(smallest_bounded, 2, 1)
    expected_smallest = {(1, 0): 2, (3, 0): 2, (3, 2): 1, (5, 0): 1}
    assert actual_decomposition(smallest_bounded) == expected_smallest
    log("bounded two-factor control K_2 S_2^2", sorted(expected_smallest.items()))

    bad = K_table(3)
    bad = multiply_factor(bad, 2, 1)
    bad = multiply_factor(bad, 2, 1)
    bad_decomp = actual_decomposition(bad)
    expected_bad = {
        (0, 0): 2, (1, 1): -2, (2, 0): 4, (2, 2): 3,
        (3, 3): -1, (4, 0): 2, (4, 2): 2, (6, 0): 1}
    assert bad_decomp == expected_bad, bad_decomp
    log("exact block-criterion obstruction K_3 S_2^2", sorted(bad_decomp.items()))
    cases, floor = verify_prefix_class()
    log("PASS two-minus prefix checks", cases, "cases; minimum prefix", floor)
    log("PASS FM-STR9e verifier")

main()
