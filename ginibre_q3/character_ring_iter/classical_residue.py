"""Compute L = Q_3(11) - 4*Q_3(9) for each small classical group in the
m_3=1 non-stable residue.

Uses the existing character-ring iteration (imported from local code).
Groups: B_2..B_8, C_2..C_7, D_5..D_8.
"""
import sys
sys.path.insert(0, '/home/sky/yang4/ginibre_q3/character_ring_iter')
import lie_data

# Generate Cartan matrices for B_n, C_n, D_n
def cartan_Bn(n):
    # rank n. simple roots α_1..α_{n-1} long (|α|²=2), α_n short (|α|²=1).
    # Cartan A[i][j] = 2(α_i, α_j)/(α_j, α_j)
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    for i in range(n-1):
        if i < n-2:
            A[i][i+1] = -1
            A[i+1][i] = -1
        else:  # i = n-2: connection between α_{n-2} (long) and α_{n-1} (short)
            A[i][i+1] = -2  # 2(α_i, α_{i+1})/(α_{i+1}, α_{i+1}) = 2·(-1)/1 = -2
            A[i+1][i] = -1  # 2(α_{i+1}, α_i)/(α_i, α_i) = 2·(-1)/2 = -1
    return A

def cartan_Cn(n):
    # rank n. simple α_1..α_{n-1} short (|α|²=1), α_n long (|α|²=2).
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    for i in range(n-1):
        if i < n-2:
            A[i][i+1] = -1
            A[i+1][i] = -1
        else:
            A[i][i+1] = -1  # short to long: 2(α_i, α_{i+1})/(α_{i+1}, α_{i+1}) = 2·(-1)/2 = -1
            A[i+1][i] = -2  # long to short: 2(α_{i+1}, α_i)/(α_i, α_i) = 2·(-1)/1 = -2
    return A

def cartan_Dn(n):
    # rank n. simply-laced. Bourbaki ordering: α_1, α_2 both adjacent to α_3 (0-indexed: 0, 1 both adjacent to 2); α_3..α_n chain.
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    # α_0 adjacent to α_2
    A[0][2] = -1; A[2][0] = -1
    # α_1 adjacent to α_2
    A[1][2] = -1; A[2][1] = -1
    # α_2, α_3, ..., α_{n-1} form a chain
    for i in range(2, n-1):
        A[i][i+1] = -1
        A[i+1][i] = -1
    return A

# Add to lie_data
for n in [2, 3, 4, 5, 6, 7, 8]:
    lie_data.CARTAN_MATRICES[f'B{n}'] = cartan_Bn(n)
for n in [2, 3, 4, 5, 6, 7]:
    lie_data.CARTAN_MATRICES[f'C{n}'] = cartan_Cn(n)
for n in [4, 5, 6, 7, 8]:
    lie_data.CARTAN_MATRICES[f'D{n}'] = cartan_Dn(n)

# Use character_ring iteration directly (Python version for small groups — fast enough for n ≤ 13)
from character_ring import compute_moments as iterate_moments

for G in ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8',
          'C2', 'C3', 'C4', 'C5', 'C6', 'C7',
          'D4', 'D5', 'D6', 'D7', 'D8']:
    try:
        moms = iterate_moments(G, 13)
        m = [moms[k] for k in range(14)]
        # Compute Q_3(n) = 2 sum C(n,k) (m_{k+2} m_{n-k} - m_{k+1} m_{n-k+1})
        from math import comb
        def Q3(n):
            return 2 * sum(comb(n, k) * (m[k+2]*m[n-k] - m[k+1]*m[n-k+1]) for k in range(n+1))
        L = Q3(11) - 4 * Q3(9)
        print(f'{G}: m=(m_0..m_6) = {m[:7]};  L = Q_3(11) - 4Q_3(9) = {L} {"✓" if L > 0 else "✗"}')
    except Exception as e:
        print(f'{G}: ERROR {e}')
