#!/usr/bin/env python3
"""E_n = E_{Sp(2n)}[prod_i S^(n)_{p_i}], S^(n)_p = sum_j chi_p(theta_j), via the Heine/Andreief
determinant E_n[prod_j F(theta_j)] = det[<F chi_j chi_k>]_{j,k<n} applied to F=exp(sum s_p chi_p),
extracting the multilinear coefficient.  Tests E_n >= 0 and monotonicity in n.  usage: probe_su2_fm3_sp2n_mono.py L K NMAX"""
import sys, itertools
from fractions import Fraction
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
L,K,NMAX=map(int,sys.argv[1:4])
def cg(a,p): return range(abs(a-p),a+p+1,2)
def mult_all(word):
    """dict a -> multiplicity of V_a in tensor of labels in word"""
    d={0:1}
    for p in word:
        nd={}
        for a,v in d.items():
            for c in cg(a,p): nd[c]=nd.get(c,0)+v
        d=nd
    return d
def En(word,n):
    # multilinear coefficient of det[G_jk], G_jk=<F chi_j chi_k>, F=prod (1+s_i chi_{p_i}) multilinear
    # E_n = sum over ordered set partitions of word into n blocks (one per row) of det-expansion:
    # det = sum_{sigma} sgn prod_j G_{j,sigma(j)}; multilinear coeff = sum over assignments of factors to rows.
    k=len(word); tot=0
    for assign in itertools.product(range(n),repeat=k):
        blocks=[[word[i] for i in range(k) if assign[i]==j] for j in range(n)]
        ms=[mult_all(b) for b in blocks]
        # G_{j,l} restricted to block j: <V_block chi_j chi_l> = sum_c m_c * N(c, j, l)
        def G(j,l,m):
            return sum(v for c,v in m.items() if (abs(j-l)<=c<=j+l and (c+j+l)%2==0))
        # det with row j using block j
        M=[[G(j,l,ms[j]) for l in range(n)] for j in range(n)]
        # determinant (small n)
        def det(M):
            if len(M)==1: return M[0][0]
            return sum((-1)**c*M[0][c]*det([row[:c]+row[c+1:] for row in M[1:]]) for c in range(len(M)))
        tot+=det(M)
    return tot
viol_pos=0; viol_mono=0; tested=0
for k in range(1,K+1):
    for word in itertools.combinations_with_replacement(range(1,L+1),k):
        vals=[En(list(word),n) for n in range(1,NMAX+1)]
        tested+=1
        if any(v<0 for v in vals): viol_pos+=1; print("NEG",word,vals)
        if any(vals[i+1]<vals[i] for i in range(len(vals)-1)):
            viol_mono+=1
            if viol_mono<=15: print("nonmonotone",word,vals)
print("words=%d negatives=%d nonmonotone=%d"%(tested,viol_pos,viol_mono))
