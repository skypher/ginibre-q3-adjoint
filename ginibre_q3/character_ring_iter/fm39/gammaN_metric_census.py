# gamma <= N branch: phi = delta - U^V (FM-SEC45 notation).  Test the single-metric bound directly with the exact delta:
#   exists t in (0, 4V):  16 t (4V - t) delta^2 >= F_U(t) F_V(t)   (then |U^V| <= delta, so phi >= 0),
# F_Y(t) = alpha_Y t + beta_Y, alpha_Y = 4(Y1^2 - Y2^2), beta_Y = 4(d Y1 - sigma Y2)^2 (FM-MECH26's Q_t, det = (4V-t)/t).
# Exact: feasibility of -K t^2 + L t - Z >= 0 on (0, 4V), K = 16 delta^2 + aU aV, L = 64 V delta^2 - aU bV - bU aV, Z = bU bV.
import sys, time
from fractions import Fraction
exec(open('gammaN_cert_lib.py').read())
from collections import Counter
def feasible(K,L,Z,T):
    f=lambda t: -K*t*t+L*t-Z
    cands=[Fraction(1,10**6),Fraction(T)-Fraction(1,10**6)]
    if K!=0:
        tv=Fraction(L,2*K)
        if 0<tv<T: cands.append(tv)
    return any(f(t)>=0 for t in cands)
EMAX,AMAX=int(sys.argv[1]),int(sys.argv[2]); AMIN=int(sys.argv[3]) if len(sys.argv)>3 else 0
st=Counter(); ex=[]; t0=time.time(); last=t0
EMIN=int(sys.argv[4]) if len(sys.argv)>4 else 3
for e in range(EMIN,EMAX+1,2):
    for a in range(AMIN,AMAX+1):
        N=a+e; V0=(a+1)*(e+1); c=c_row(a,e); d0,sigma=a-e,N+2
        for C in range(3,N+1):
            for s_ in range(0,(N-C)//2+1):
                rest=N-C-2*s_
                for q in range(0,rest//2+1):
                    p=rest-q
                    alpha,tau=p+s_,N-s_+1
                    U=add(g(c,alpha),g(c,q+s_)); Vv=add(g(c,tau),g(c,s_-C-1)); X=wedge(U,Vv)
                    delta=P(c,alpha,C)-P(c,s_-C-1,C)
                    st['words']+=1
                    if delta-X<0: st['PHI<0']+=1
                    if X==0 or (delta>=abs(X) and False): pass
                    au=4*(U[0]**2-U[1]**2); bu=4*(d0*U[0]-sigma*U[1])**2
                    av=4*(Vv[0]**2-Vv[1]**2); bv=4*(d0*Vv[0]-sigma*Vv[1])**2
                    K=16*delta**2+au*av; L=64*V0*delta**2-au*bv-bu*av; Z=bu*bv
                    if delta>=0 and feasible(K,L,Z,4*V0): st['metric']+=1
                    elif X<=0 and delta>=0: st['X<=0 (trivial)']+=1
                    else:
                        st['NOT metric']+=1
                        if len(ex)<12: ex.append((e,a,C,p,q,s_,delta,X))
        if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat e=%d a=%d'%(e,a),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed %.1f'%(time.time()-t0)); print('not-metric examples (e,a,C,p,q,s,delta,X):',ex)
