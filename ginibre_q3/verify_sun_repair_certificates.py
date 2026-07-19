#!/usr/bin/env python3
"""Legacy Python cross-check for the SU_N_REPAIR endpoint certificates.

The active load-bearing replay is the parallel GMP implementation
``character_ring_iter/verify_sun_repair_certificates_gmp.cpp``.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from time import time


def partitions_bounded_length(n: int, max_part: int, max_len: int):
    if n == 0:
        yield ()
        return
    if max_len == 0:
        return
    for first in range(min(max_part, n), 0, -1):
        for rest in partitions_bounded_length(n - first, first, max_len - 1):
            yield (first,) + rest


def hook_dim(lam: tuple[int, ...], facts: list[int]) -> int:
    ell = len(lam)
    n = sum(lam)
    num = facts[n]
    den = 1
    for i, li in enumerate(lam):
        den *= facts[li + ell - i - 1]
    prod = 1
    for i in range(ell):
        for j in range(i + 1, ell):
            prod *= lam[i] - lam[j] + j - i
    return num * prod // den


def trace_moments(N: int, max_r: int, progress: bool = False) -> list[Fraction]:
    facts = [factorial(i) for i in range(max_r + 1)]
    out: list[Fraction] = []
    start = time()
    for r in range(max_r + 1):
        total = 0
        count = 0
        for lam in partitions_bounded_length(r, r, N):
            dim = hook_dim(lam, facts)
            total += dim * dim
            count += 1
        out.append(Fraction(total))
        if progress and (r % 10 == 0 or r == max_r):
            print(
                f"moments N={N} r={r}/{max_r} partitions={count} elapsed={time()-start:.1f}s",
                flush=True,
            )
    return out


def adjoint_moments(E: list[Fraction], max_r: int) -> list[Fraction]:
    return [
        sum(Fraction(comb(r, a)) * ((-1) ** (r - a)) * E[a] for a in range(r + 1))
        for r in range(max_r + 1)
    ]


def shifted_moment(E: list[Fraction], n: int, center: Fraction) -> Fraction:
    return sum(Fraction(comb(n, t)) * (-center) ** (n - t) * E[t] for t in range(n + 1))


def prop64_prop65() -> None:
    delta_b_64 = Fraction(244495424, 25 * 1065343**2)
    delta_c_64 = Fraction(471756427, 178410035553750000)
    value_64 = delta_b_64 * delta_c_64 * Fraction(3, 2) ** 89
    assert value_64 > 72

    eta_b_65 = Fraction(122923229137548665407, 65536 * 1467426013**2)
    eta_c_65 = Fraction(88255484124721, 992717442773183102976)
    value_65 = eta_b_65 * eta_c_65 * Fraction(3, 2) ** 71
    assert value_65 > 200
    print("Propositions 64-65 endpoint constants: ok", flush=True)


def prop66() -> None:
    rows = {
        24: range(57, 70, 2),
        25: range(59, 70, 2),
        26: range(61, 70, 2),
        27: range(65, 70, 2),
    }
    for N, ks in rows.items():
        vals: list[tuple[Fraction, int]] = []
        for k in ks:
            K = k + 2
            tau = Fraction(9 * comb(K, N + 1), factorial(N + 1))
            bracket = Fraction(7, 18) - tau * (4 * k + 10) - tau * tau * (2 * k + 5)
            vals.append((bracket, k))
        min_bracket, min_k = min(vals)
        assert min_bracket > Fraction(3, 10)
        print(f"Proposition 66 row N={N}: min at k={min_k}, ok", flush=True)


def prop67() -> None:
    coeff = [959324176457, 276017837985, -58083597219, -11664458272, 1230335979]
    y_lo = Fraction(29, 10)
    y_hi = Fraction(57, 10)
    center = Fraction(43, 10)
    radius = Fraction(7, 5)
    S = Fraction(935204207737834, 625)

    band_target = Fraction(44525786465344542780431067563007, 10_000_000_000)
    delta_b = band_target / (Fraction(49, 25) * S * S)
    assert delta_b == Fraction(
        44525786465344542780431067563007,
        43884276324717505323728973379833856,
    )

    min_excess: tuple[Fraction, int] | None = None
    for N in range(6, 10):
        E = trace_moments(N, 10)
        value = Fraction(0)
        for i, ci in enumerate(coeff):
            for j, cj in enumerate(coeff):
                n = i + j
                h_moment = radius * radius * shifted_moment(E, n, center) - shifted_moment(
                    E, n + 2, center
                )
                value += Fraction(ci) * Fraction(cj) * h_moment
        excess = value - band_target
        assert excess >= 0
        if min_excess is None or excess < min_excess[0]:
            min_excess = (excess, N)
    assert min_excess == (Fraction(1513726621221888441), 9)
    print("Proposition 67 band certificate: ok", flush=True)

    delta_c = Fraction(471756427, 178410035553750000)
    tail_min: tuple[Fraction, int] | None = None
    for N in range(6, 15):
        E = trace_moments(N, 30)
        tail = (1 - Fraction(6**15, E[15])) ** 2 * E[15] * E[15] / E[30]
        excess = tail - delta_c
        assert excess >= 0
        if tail_min is None or excess < tail_min[0]:
            tail_min = (excess, N)
    assert tail_min == (
        Fraction(
            2602599137945210797524667371658649,
            56336839730247629366109305276233265935886250000,
        ),
        14,
    )

    endpoint = delta_b * delta_c * Fraction(19, 10) ** 52 * Fraction(3, 10) ** 2
    assert endpoint > 72
    print("Proposition 67 tail and endpoint certificates: ok", flush=True)


def prop68() -> None:
    expected = {
        6: (
            39,
            int(
                "64535719854392029544773326762011774756925778563684823690382903928732901693552607059275"
            ),
        ),
        7: (
            41,
            int(
                "81983250005769227285498956748387927182179586638327455931621250522805530412550262461119689354176"
            ),
        ),
        8: (
            41,
            int(
                "13776730806844857067459723475572806115606561762537749689970379383403774177954665948696550730822960"
            ),
        ),
        9: (
            43,
            int(
                "516741713808588211599059495979604717390970752707423594032209037987274754148154730954160916056953963450368"
            ),
        ),
        10: (
            43,
            int(
                "3618232844122202584832213951704951536993240226634987294481850293593859749135402230512808610277609410270817"
            ),
        ),
        11: (
            45,
            int(
                "33399774449175848055011934469656070323910964218709116241112806497510841469515062121126196910235665905718776943228"
            ),
        ),
        12: (
            45,
            int(
                "60567314448754365018106205779950325082006402099956453971960484870790769677108940276606446856248401017048182106706"
            ),
        ),
        13: (
            47,
            int(
                "363239305623168787675421494307053710143238482046635190860778910231693491256558565289022469174629997636148032811272314736"
            ),
        ),
        14: (
            47,
            int(
                "412017453922855017028903546300156146196642966415729406173097639746151298064337858733695094575907116122702045750587973615"
            ),
        ),
        15: (
            49,
            int(
                "2454353316530929119018984347889491575535943771589993590419065745590973632964833686716148957357041276645873773505319638231954239"
            ),
        ),
        16: (
            49,
            int(
                "2493119549736314583261850567832642447348332533596249440198341404887202602996116002199771825733723134557950156991608933086839975"
            ),
        ),
        17: (
            51,
            int(
                "16910085949149500902700908171041557824229854633053201817574276758325371717389977268687410623596430698509784622587525716043509684517756"
            ),
        ),
        18: (
            51,
            int(
                "16928441762086713834411725710317433241255361224654439216646167863265666371228166615581072054953650607845481890940623126734747743109775"
            ),
        ),
    }
    start = time()
    for N, (expected_k, expected_delta) in expected.items():
        first_k = N + 33
        if first_k % 2 == 0:
            first_k += 1
        ks = list(range(first_k, 53, 2))
        E = trace_moments(N, max(ks) + 2, progress=True)
        m = adjoint_moments(E, max(ks) + 2)
        got_delta, got_k = min((m[k + 2] * m[k] - m[k + 1] * m[k + 1], k) for k in ks)
        assert (got_k, int(got_delta)) == (expected_k, expected_delta)
        print(f"Proposition 68 row N={N}: ok, elapsed={time()-start:.1f}s", flush=True)


def main() -> None:
    prop64_prop65()
    prop66()
    prop67()
    prop68()
    print("All SU(N) repair certificates verified.", flush=True)


if __name__ == "__main__":
    main()
