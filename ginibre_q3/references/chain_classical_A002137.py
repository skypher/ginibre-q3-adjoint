"""Chain Inequality verification using OEIS A002137 stable SO(oo)/Sp(oo)
moments. These apply to SO(N), Sp(N) for N large enough that
Rains/Diaconis-Shahshahani moment stabilization has kicked in.
"""
from math import comb
from pathlib import Path

REFS = Path(__file__).resolve().parent

def read_bfile(path):
    moments = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            n_str, a_str = line.split()
            moments.append(int(a_str))
    return moments

def Q3(n, m):
    return 2 * sum(comb(n, k) * (m[k+2] * m[n-k] - m[k+1] * m[n-k+1])
                   for k in range(n+1))

def main():
    mom = read_bfile(REFS / "oeis_A002137_stable.txt")
    m_max = len(mom) - 1
    m_bar_max = (m_max - 5) // 2
    print(f"Stable SO(oo)/Sp(oo) A002137: {m_max+1} moments (m_0..m_{m_max})")
    print(f"Chain m = 0..{m_bar_max}")
    print()
    # Summarize every 10 m values
    for mb in range(0, m_bar_max + 1):
        q1 = Q3(2*mb + 1, mom)
        q3 = Q3(2*mb + 3, mom)
        diff = q3 - 4*q1
        status = 'OK' if diff >= 0 else 'FAIL'
        if mb <= 5 or mb % 5 == 0 or mb == m_bar_max:
            print(f"  Chain m={mb:3d}: diff = {diff} [{status}]")
    print()
    # Also Riordan for SO(3) = A_1
    riordan = read_bfile(f"{REFS}/oeis_A005043_riordan.txt") if False else None

main()
