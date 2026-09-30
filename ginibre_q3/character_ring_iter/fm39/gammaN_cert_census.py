# Coverage census of FM-SEC45's gamma <= N certificate: phi = delta - U^V, delta = sum_d (E + eps W) >= Lambda (OL or M_i/M_j
# step bounds), |U^V| <= Lambda via one Q_t metric.  Status per word: 'cert', 'step' (an adverse step has no admissible
# M_i/M_j bound), 'KLZ' (Lambda too small for the wedge).  Also records phi >= 0 (exact).
import sys, time
exec(open('gammaN_cert_lib.py').read())
from collections import Counter
def status(e,a,C,p,q,s):
    N=a+e; A,B0=C+p,C+q; V0=(a+1)*(e+1); c=c_row(a,e)
    alpha,tau=p+s,N-s+1
    U=add(g(c,alpha),g(c,q+s)); V=add(g(c,tau),g(c,s-C-1)); X=wedge(U,V)
    delta=P(c,alpha,C)-P(c,s-C-1,C)
    Lam=0
    for d in range(abs(A-B0),A+B0+1,2):
        j0=(N+d-C)//2; i=j0+C+1
        j,eps=(N-j0,1) if d<C else (j0,-1)
        E=D(c,j)-D(c,i); W=at(c,j)*B(c,i)-B(c,j)*at(c,i)
        if eps*W>=0: Lam+=E; continue
        x,y=2*j-N,2*i-N; cands=[]
        if 0<y*y<4*V0:
            rhs=D(c,i)*(4*y*y*D(c,j)-(y*y-x*x)*Ac(c,j)**2)
            if rhs>=0:
                R=Fraction(rhs,4*V0-y*y)
                if E*E>=R: cands.append(R)
        if 0<x*x<4*V0:
            rhs=D(c,j)*(4*x*x*D(c,i)+(y*y-x*x)*Ac(c,i)**2)
            if rhs>=0:
                R=Fraction(rhs,4*V0-x*x)
                if E*E>=R: cands.append(R)
        if not cands: return 'step',delta-X
        Lam+=E-ceil_sqrt(min(cands))
    d0,sigma=a-e,N+2
    au=4*(U[0]**2-U[1]**2); bu=4*(d0*U[0]-sigma*U[1])**2
    av=4*(V[0]**2-V[1]**2); bv=4*(d0*V[0]-sigma*V[1])**2
    K=16*Lam**2+au*av; L=64*V0*Lam**2-au*bv-bu*av; Z=bu*bv
    ok=K>0 and 0<L<8*V0*K and L*L>=4*K*Z
    return ('cert' if ok else 'KLZ'),delta-X
EMAX,AMAX=int(sys.argv[1]),int(sys.argv[2])
st=Counter(); ex={'step':[],'KLZ':[]}; t0=time.time(); last=t0
for e in range(1,EMAX+1,2):
    for a in range(0,AMAX+1):
        N=a+e
        for C in range(3,N+1):
            for s_ in range(0,(N-C)//2+1):
                rest=N-C-2*s_
                for q in range(0,rest//2+1):
                    p=rest-q
                    stt,phi=status(e,a,C,p,q,s_)
                    st[stt]+=1; st['words']+=1
                    if phi<0: st['PHI<0']+=1
                    if stt in ex and len(ex[stt])<6: ex[stt].append((e,a,C,p,q,s_))
        if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat e=%d a=%d'%(e,a),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed %.1f'%(time.time()-t0)); print('examples',ex)
