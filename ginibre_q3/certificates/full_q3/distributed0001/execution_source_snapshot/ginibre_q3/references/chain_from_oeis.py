"""Compute Chain Inequality Q_3(2m+3) - 4*Q_3(2m+1) for each exceptional group
using exact moments from OEIS b-files.

OEIS identifications (all by Bruce Westbury, contributed 2010-2013):
  G_2: A227292 (moment sequence for G_2 adjoint)
  F_4: A179685
  E_6: A179684
  E_7: A179683
  E_8: A179663

All five are "dimension of invariant tensors in adj^⊗n" for the respective
exceptional Lie algebra — exactly our m_k = mult(triv, adj^⊗k).
"""
from math import comb
from pathlib import Path

REFS = Path(__file__).resolve().parent

GROUPS = {
    'G_2': 'oeis_A227292.txt',
    'F_4': 'oeis_A179685.txt',
    'E_6': 'oeis_A179684.txt',
    'E_7': 'oeis_A179683.txt',
    'E_8': 'oeis_A179663.txt',
}

def read_bfile(path):
    """Read OEIS b-file: lines of 'n a(n)'."""
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

def chain_for(group, moments):
    m_max = len(moments) - 1
    # Chain m uses Q_3(2m+1) and Q_3(2m+3); Q_3(n) needs m_{n+2}.
    # So Chain m feasible iff 2m+5 ≤ m_max.
    m_bar_max = (m_max - 5) // 2
    print(f"\n{group}: {m_max+1} moments available (m_0..m_{m_max}) → "
          f"Chain m=0..{m_bar_max}")
    print(f"  moments: {moments[:6]}... (last: m_{m_max}={moments[-1]})")
    for mb in range(0, m_bar_max + 1):
        q1 = Q3(2*mb + 1, moments)
        q3 = Q3(2*mb + 3, moments)
        diff = q3 - 4*q1
        status = 'OK' if diff >= 0 else 'FAIL'
        print(f"  Chain m={mb:2d}: Q_3({2*mb+3}) - 4 Q_3({2*mb+1}) = {diff} [{status}]")

def main():
    for group, fname in GROUPS.items():
        moments = read_bfile(REFS / fname)
        chain_for(group, moments)

if __name__ == "__main__":
    main()
