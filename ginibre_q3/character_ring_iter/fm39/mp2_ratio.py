# Sp(4) multiplicities of E_kappa = (x)_i Sym^(kappa_i) W (restricted to SU(2)^2: Sym^k W = H_k):
#   m_00 = phi_1(E), m_11 = phi_1(E (xy+1)), m_20 = phi_1(E (x^2+y^2+xy-2))   (phi_1 = (1/2) E[(x-y)^2 .]).
# Sector inequality: 3 m_11 <= m_20 + 5 m_00.  Report max ratio 3 m_11/(m_20 + 5 m_00), equality cases.
exec(open('mech28_eval.py').read())
from fractions import Fraction
def partitions(n, maxpart=None):
    if maxpart is None: maxpart=n
    if n==0: yield (); return
    for p in range(min(n,maxpart),0,-1):
        for rest in partitions(n-p,p): yield (p,)+rest
xy1={(1,1):1,(0,0):1}; c20={(2,0):1,(0,2):1,(1,1):1,(0,0):-2}
best=(Fraction(0),None); eq=[]; zeros=0; n=0
for size in range(1,13):
    for kap in partitions(size):
        E=word(kap,())
        m00=phi(1,E); m11=phi(1,mul(E,xy1)); m20=phi(1,mul(E,c20))
        assert m00>=0 and m11>=0 and m20>=0
        n+=1; den=m20+5*m00
        assert 3*m11<=den, kap
        if den==0: zeros+=1; continue
        r=Fraction(3*m11)/den
        if r>best[0]: best=(r,kap,(m00,m11,m20))
        if r==1 and len(eq)<12: eq.append((kap,(m00,m11,m20)))
print('partitions',n,'with m20+5m00=0:',zeros); print('max 3m11/(m20+5m00) =',best[0],'=',float(best[0]),'at',best[1],best[2]); print('equality cases:',eq)
