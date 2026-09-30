# Quotient slot injection (QGSI) test for EM2.
# EM2 at kappa = (kappa'', a, b) (a>=b the two smallest parts, last two slots) is T1 for
# Mbar = h_kappa'' (x) D_(a,b), realised as ker(Delta) where Delta = omega-contraction of the last two slots.
# Slot maps A_i, B_i for i < s-2 commute with Delta, so they map Hom(U,Mbar) -> Hom(1,Mbar), Hom(ad,Mbar).
# Test: generic T-bar = (5 A-rows, 1 B-row) on 3 copies of Hom(U,Mbar) has rank 3*mbar_U.
import os, sys
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'slotHWg.py')).read().split("if __name__=='__main__':")[0]
exec(src)
def delta_op(W,Wp,w):
    # Delta: M_kappa (weight w) -> M_kappa' (weight w), kappa' = kappa - e_{s-2} - e_{s-1} (0-based last two)
    s=W.s; src_=W.I(w); dst=Wp.I(w); st=W.states(w); n=len(src_)
    ex1=W.bases[s-2][st[s-2]]; ex2=W.bases[s-1][st[s-1]]
    idx1={tuple(b):i for i,b in enumerate(Wp.bases[s-2].tolist())}
    idx2={tuple(b):i for i,b in enumerate(Wp.bases[s-1].tolist())}
    # destination flat index: same for other slots, replaced for last two
    other=np.zeros(n,dtype=np.int64)
    for i in range(s-2): other+=st[i]*Wp.strides[i]
    R=[];C=[];V=[]
    for (u,v_),sg in (((0,2),1),((2,0),-1),((1,3),1),((3,1),-1)):
        ok=(ex1[:,u]>0)&(ex2[:,v_]>0)
        cols=np.flatnonzero(ok)
        for c in cols.tolist():
            b1=ex1[c].copy(); b1[u]-=1; b2=ex2[c].copy(); b2[v_]-=1
            tgt=other[c]+idx1[tuple(b1.tolist())]*Wp.strides[s-2]+idx2[tuple(b2.tolist())]*Wp.strides[s-1]
            r=np.searchsorted(dst,tgt); assert dst[r]==tgt
            R.append(r); C.append(c); V.append((sg*int(ex1[c,u])*int(ex2[c,v_]))%P)
    if not R: return sps.csr_matrix((len(dst),n),dtype=np.int64)
    M=sps.csr_matrix((V,(R,C)),shape=(len(dst),n),dtype=np.int64); M.data%=P; M.eliminate_zeros(); return M
def runq(kap,seed=5):
    t0=time.time(); STATE['word']=kap; s=len(kap)
    kp=kap[:-2]+(kap[-2]-1,kap[-1]-1)
    m=char_mults(kap,[(0,0),(1,1),(2,0)]); mp=char_mults(kp,[(0,0),(1,1),(2,0)])
    m1,mU,mad=[x-y for x,y in zip(m,mp)]
    log(f'{kap}: quotient mults (m0,mU,mad)=({m1},{mU},{mad}); 3mU={3*mU}, 5m0+mad={5*m1+mad}, EM2 margin={5*m1+mad-3*mU}')
    if mU==0: log('  trivially injective'); STATE['word']=None; return True,dict(m1=m1,mU=mU,mad=mad,rank=0)
    W=Word(kap); Wp=Word(kp); n11=len(W.I((1,1)))
    stage('root kernels')
    Ks={r_:ker_root(W,r_) for r_ in (1,2)}
    root=min((1,2),key=lambda r_:Ks[r_].shape[1]); K2=Ks[root]; n2=K2.shape[1]; del Ks
    other,wo=((ea1,(2,0)) if root==2 else (ea2,(1,3)))
    R1=W.op(other,(1,1),wo); A1=(R1@K2).tocsr(); A1.data%=P; A1.eliminate_zeros()
    A1=A1[np.flatnonzero(np.diff(A1.indptr))]
    A=fill_nmod(A1); del A1
    rk,N=nullspace_inplace(A); del A
    if N.shape[0]!=m[1]: log(f'  full HW nullity {N.shape[0]} != m_U {m[1]}'); return False,{}
    # restrict to ker Delta
    D=delta_op(W,Wp,(1,1))
    DK=(D@K2).tocsr(); DK.data%=P
    Z=np.asarray(DK.todense(),dtype=np.int64) if D.shape[0] else np.zeros((0,n2),dtype=np.int64)
    Z=mulmod(Z%P,(N.T)%P) if Z.shape[0] else Z
    if Z.shape[0]:
        Az=flint.nmod_mat(Z.shape[0],Z.shape[1],(Z%P).ravel().tolist(),P)
        Cn=Az.nullspace(); Cm,kz=Cn
        C=np.array([[int(Cm[i,j]) for j in range(kz)] for i in range(Z.shape[1])],dtype=np.int64).reshape(Z.shape[1],kz)
    else:
        C=np.eye(N.shape[0],dtype=np.int64); kz=N.shape[0]
    if kz!=mU: log(f'  harmonic HW dim {kz} != mbar_U {mU}  MISMATCH'); return False,{}
    Nq=mulmod(C.T%P,N%P)   # mbar_U x n2
    rng=np.random.default_rng(seed)
    slots=list(range(s-2)) if os.environ.get('SLOTS','inner')=='inner' else list(range(s))
    a=rng.integers(1,P,size=(5,3,s)); b=rng.integers(1,P,size=(3,s))
    n00=len(W.I((0,0))); n20=len(W.I((2,0)))
    dA=min(n00,m[0]+8); dB=min(n20,m[2]+8)
    Lops=[W.op(L,u,v) for L,(u,v) in zip((fa2,fa1,fa1,fa2),zip(Ywt[:-1],Ywt[1:]))]
    SA=[[W.op(Yd[k],Ywt[k],(0,0),slots=[i]) for k in range(5)] for i in range(s)]
    SB=[(W.op(Y[0],(1,-1),(2,0),slots=[i]),W.op(Y[1],(1,1),(2,0),slots=[i])) for i in range(s)]
    def proj(d,mm,nnzc=24):
        r_=rng.integers(0,d,size=mm*nnzc); c_=np.repeat(np.arange(mm),nnzc); v_=rng.integers(1,P,size=mm*nnzc)
        M_=sps.csr_matrix((v_,(r_,c_)),shape=(d,mm),dtype=np.int64); M_.data%=P; return M_
    PA=proj(dA,n00); PB=proj(dB,n20)
    for Q_ in [K2]+Lops:
        if Q_.nnz: assert int(np.diff(Q_.tocsr().indptr).max())*(P-1)**2 < 2**62
    F=[(K2@Nq.T)%P]
    for L in Lops: F.append((L@F[-1])%P)
    TA=np.zeros((3,5,dA,mU),dtype=np.int64); TB=np.zeros((3,dB,mU),dtype=np.int64)
    def red(Msp):
        Msp=Msp.tocsr(); Msp.data%=P; Msp.eliminate_zeros(); return Msp
    for i in slots:
        # every sparse factor is reduced mod P before the dense product (entries < P, row nnz small: no int64 overflow)
        QAi=[red(PA@SA[i][k]) for k in range(5)]; QB0=red(PB@SB[i][0]); QB1=red(PB@SB[i][1])
        for Q_ in QAi+[QB0,QB1]:
            if Q_.nnz: assert int(np.diff(Q_.indptr).max())*(P-1)**2 < 2**62, 'row nnz too large for exact int64 product'
        x=(QAi[0]@F[0])%P
        for k in range(1,5): x=(x+(QAi[k]@F[k])%P)%P
        y=((QB0@F[1])%P-(QB1@F[0])%P)%P
        for c in range(3):
            for r in range(5): TA[c,r]=(TA[c,r]+int(a[r,c,i])*x)%P
            TB[c]=(TB[c]+int(b[c,i])*y)%P
    Tm=np.concatenate([np.concatenate([TA[c,r] for r in range(5)]+[TB[c]],axis=0).T for c in range(3)],axis=0)
    rkT=rank_dense(Tm)
    if os.environ.get('CHECKH','1')=='1' and len(W.I((0,0)))<=200000:
        D0=delta_op(W,Wp,(0,0)); worst=0
        for i in slots:
            x=(red(SA[i][0])@F[0])%P
            for k in range(1,5): x=(x+(red(SA[i][k])@F[k])%P)%P
            z=(red(D0)@x)%P; worst=max(worst,int(np.count_nonzero(z)))
        log(f'  harmonic check: nonzero entries of Delta(A_i images) = {worst}')
    # also single-copy ranks: rank of A-part alone and B-part alone
    ok=rkT==3*mU
    log(f'  rank(Tbar)={rkT} / 3mbarU={3*mU} [slots={"inner" if len(slots)==s-2 else "all"}]: {"INJECTIVE" if ok else "NOT injective"} ({time.time()-t0:.1f}s)')
    STATE['word']=None
    return ok,dict(m1=m1,mU=mU,mad=mad,rank=rkT)
if __name__=='__main__':
    for arg in sys.argv[1:]:
        runq(tuple(int(x) for x in arg.split(',')))
