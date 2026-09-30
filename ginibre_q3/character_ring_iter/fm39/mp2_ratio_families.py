# Trend of 3 m_11/(m_20 + 5 m_00) for E_kappa along families (Sp(4) multiplicities via SU(2)^2 characters).
exec(open('mech28_eval.py').read())
from fractions import Fraction
xy1={(1,1):1,(0,0):1}; c20={(2,0):1,(0,2):1,(1,1):1,(0,0):-2}
def ratio(kap):
    E=word(kap,()); m00=phi(1,E); m11=phi(1,mul(E,xy1)); m20=phi(1,mul(E,c20))
    den=m20+5*m00
    return (float(Fraction(3*m11)/den) if den else None), (m00,m11,m20)
fams={'(1^n)':lambda k:(1,)*(2*k),'(2^n)':lambda k:(2,)*k,'(3^k,1^k)':lambda k:(3,)*k+(1,)*k,'(3^n)':lambda k:(3,)*(2*k),'(2,1^n)':lambda k:(2,)+(1,)*(2*k),'(4^k,2^k)':lambda k:(4,)*k+(2,)*k}
for name,f in fams.items():
    out=[]
    for k in range(1,9):
        kap=f(k)
        if sum(kap)>36: break
        r,_=ratio(kap); out.append((k, None if r is None else round(r,4)))
    print(name, out, flush=True)
