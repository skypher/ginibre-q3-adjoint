# Memory-lean highest-weight slot-map test for T1.
# Same test as slotHWf.py (generic slot-map T on Hom(U,M)^3, rank vs 3m_U), but:
#  - ker e_a2 on M_(1,1) is built block by block (e_a2 = x1 d/dx3 preserves each slot's (e0,e2)),
#    so the dense elimination only sees ker e_a2 (about a third of M_(1,1));
#  - flint matrices are filled from sparse entries and reduced in place (no dense numpy / list copies);
#  - m_1, m_U, m_ad come from the Weyl character (m_U is cross-checked against the nullity);
#  - A_i, B_i images are randomly projected to m_1+8 / m_ad+8 coordinates before the rank test
#    (a projection can only lower rank, so full rank after projection still certifies injectivity).
import numpy as np, itertools, sys, time, os, threading, resource
import scipy.sparse as sps
import flint
P=1000003
flint.ctx.threads=int(os.environ.get('FLINT_THREADS','8'))
def log(*a): print(time.strftime('[%H:%M:%S]'),*a,flush=True)
def peak_gb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6
STATE={'word':None,'stage':'idle','t0':time.time()}
def stage(s): STATE['stage']=s; STATE['ts']=time.time()
def heartbeat(period):
    while True:
        time.sleep(period)
        if STATE['word'] is not None:
            log(f'HEARTBEAT word={STATE["word"]} stage={STATE["stage"]} '
                f'stage_elapsed={time.time()-STATE.get("ts",time.time()):.0f}s peakRSS={peak_gb():.2f}GB')
def E(a,b):
    X=np.zeros((4,4),dtype=np.int64); X[a,b]=1; return X
def omega(u,w): return {(0,2):1,(1,3):1,(2,0):-1,(3,1):-1}.get((u,w),0)
def Xuv(u,v):
    X=np.zeros((4,4),dtype=np.int64)
    for w in range(4): X[v,w]+=omega(u,w); X[u,w]-=omega(v,w)
    return X
ea1=E(0,1)-E(3,2); ea2=E(1,3); fa1=E(1,0)-E(2,3); fa2=E(3,1)
Y=[Xuv(0,1)]
for L in (fa2,fa1,fa1,fa2): Y.append(L@Y[-1]-Y[-1]@L)
Ywt=[(1,1),(1,-1),(0,0),(-1,1),(-1,-1)]
G=np.array([[int(np.trace(a@b)) for b in Y] for a in Y],dtype=object)
import sympy as sp
Ginv=sp.Matrix(G).inv()
def modq(x): x=sp.Rational(x); return (int(x.p)*pow(int(x.q),P-2,P))%P
Yd=[]
for k in range(5):
    D=np.zeros((4,4),dtype=np.int64)
    for j in range(5): D=(D+modq(Ginv[k,j])*Y[j])%P
    Yd.append(D)
def symbasis(k): return [e for e in itertools.product(range(k+1),repeat=4) if sum(e)==k]
def glsym(k,X):
    B=symbasis(k); idx={b:i for i,b in enumerate(B)}; n=len(B); rows=[];cols=[];vals=[]
    for j,b in enumerate(B):
        for a in range(4):
            for c in range(4):
                if X[a,c]%P==0 or b[c]==0: continue
                nb=list(b); nb[c]-=1; nb[a]+=1
                rows.append(idx[tuple(nb)]); cols.append(j); vals.append((int(X[a,c])*b[c])%P)
    return sps.csc_matrix((vals,(rows,cols)),shape=(n,n),dtype=np.int64)
# ---------- characters ----------
WEYL=[]
for perm in ((0,1),(1,0)):
    for sg in itertools.product((1,-1),repeat=2):
        sgn=(1 if perm==(0,1) else -1)*sg[0]*sg[1]
        WEYL.append((perm,sg,sgn))
def weyl_act(w,v): perm,sg,_=w; return (sg[0]*v[perm[0]],sg[1]*v[perm[1]])
def char_mults(kap,lams):
    from collections import Counter
    C=Counter({(0,0):1})
    for k in kap:
        S=Counter()
        for e in symbasis(k): S[(e[0]-e[2],e[1]-e[3])]+=1
        D=Counter()
        for a,x in C.items():
            for b,y in S.items(): D[(a[0]+b[0],a[1]+b[1])]+=x*y
        C=D
    out=[]
    for lam in lams:
        m=0
        for w in WEYL:
            wr=weyl_act(w,(2,1)); m+=w[2]*C.get((lam[0]+2-wr[0],lam[1]+1-wr[1]),0)
        out.append(m)
    return out
# ---------- sparse operators between weight spaces ----------
class Word:
    def __init__(self,kap):
        self.kap=kap; self.s=len(kap)
        self.bases=[np.array(symbasis(k),dtype=np.int64) for k in kap]
        self.dims=[len(b) for b in self.bases]; self.dimM=int(np.prod(self.dims))
        st=[1]*self.s
        for i in range(self.s-2,-1,-1): st[i]=st[i+1]*self.dims[i+1]
        self.strides=np.array(st,dtype=np.int64)
        wx=np.zeros(1,dtype=np.int16); wy=np.zeros(1,dtype=np.int16)
        for b in self.bases:
            bx=(b[:,0]-b[:,2]).astype(np.int16); by=(b[:,1]-b[:,3]).astype(np.int16)
            wx=(wx[:,None]+bx[None,:]).ravel(); wy=(wy[:,None]+by[None,:]).ravel()
        self.wx=wx; self.wy=wy; self._I={}; self._st={}; self._G={}
    def I(self,w):
        if w not in self._I: self._I[w]=np.nonzero((self.wx==w[0])&(self.wy==w[1]))[0].astype(np.int64)
        return self._I[w]
    def states(self,w):
        if w not in self._st: self._st[w]=[a.astype(np.int64) for a in np.unravel_index(self.I(w),self.dims)]
        return self._st[w]
    def G(self,i,X):
        key=(self.kap[i],X.tobytes())
        if key not in self._G: self._G[key]=glsym(self.kap[i],X)
        return self._G[key]
    def op(self,X,wsrc,wdst,slots=None):
        src=self.I(wsrc); dst=self.I(wdst); st=self.states(wsrc); n=len(src)
        R=[];C=[];V=[]
        for i in (range(self.s) if slots is None else slots):
            Gm=self.G(i,X); b=st[i]
            cnt=np.diff(Gm.indptr)[b]; tot=int(cnt.sum())
            if tot==0: continue
            col=np.repeat(np.arange(n,dtype=np.int64),cnt)
            off=np.arange(tot,dtype=np.int64)-np.repeat(np.cumsum(cnt)-cnt,cnt)
            pos=Gm.indptr[b][col]+off
            nb=Gm.indices[pos].astype(np.int64); val=Gm.data[pos]
            tgt=src[col]+(nb-b[col])*self.strides[i]
            r=np.minimum(np.searchsorted(dst,tgt),len(dst)-1)
            assert (dst[r]==tgt).all(), 'operator leaves target weight space'
            R.append(r);C.append(col);V.append(val)
        if not R: return sps.csr_matrix((len(dst),n),dtype=np.int64)
        M=sps.csr_matrix((np.concatenate(V),(np.concatenate(R),np.concatenate(C))),shape=(len(dst),n),dtype=np.int64)
        M.data%=P; M.eliminate_zeros(); return M
# ---------- ker of one simple-root raising operator on M_(1,1), block by block ----------
# e_a2 = x1 d/dx3 keeps each slot's (e0,e2): one sl2 component per slot, (top x1, m=e1+e3), weight 1.
# e_a1 = x0 d/dx1 - x3 d/dx2 keeps each slot's (e0+e1, e2+e3): two components per slot,
#        (top x0, m=e0+e1, sign +) and (top x3, m=e2+e3, sign -), weight 0.
_loc_cache={}
def local_kernel(mt,sg,w):
    # tensor product of sl2 components Sym^{m_j}; raising t_j -> t_j+1 with coefficient sg_j*(m_j-t_j)
    key=(mt,sg,w)
    if key in _loc_cache: return _loc_cache[key]
    src=[e for e in itertools.product(*[range(m+1) for m in mt]) if sum(2*x-m for x,m in zip(e,mt))==w]
    tgt=[e for e in itertools.product(*[range(m+1) for m in mt]) if sum(2*x-m for x,m in zip(e,mt))==w+2]
    ti={e:j for j,e in enumerate(tgt)}
    if not tgt:
        K=np.eye(len(src),dtype=np.int64)
    else:
        A=flint.nmod_mat(len(tgt),len(src),P)
        for c,e in enumerate(src):
            for i,(x,m) in enumerate(zip(e,mt)):
                if m-x>0:
                    ne=list(e); ne[i]+=1; A[ti[tuple(ne)],c]=(sg[i]*(m-x))%P
        X,k=A.nullspace()
        K=np.array([[int(X[i,j]) for i in range(len(src))] for j in range(k)],dtype=np.int64).reshape(k,len(src))
    _loc_cache[key]=(len(src),K); return _loc_cache[key]
def ker_root(W,root):
    st=W.states((1,1)); n=len(W.I((1,1))); s=W.s
    ex=[W.bases[i][st[i]] for i in range(s)]          # n x 4 exponents per slot
    if root==2:
        keys=np.stack([ex[i][:,0]*(W.kap[i]+1)+ex[i][:,2] for i in range(s)],axis=1)
        top=np.stack([ex[i][:,1] for i in range(s)],axis=1)
        mts=np.stack([ex[i][:,1]+ex[i][:,3] for i in range(s)],axis=1)
        sg=(1,)*s; w=1
    else:
        keys=np.stack([ex[i][:,0]+ex[i][:,1] for i in range(s)],axis=1)
        top=np.stack([c for i in range(s) for c in (ex[i][:,0],ex[i][:,3])],axis=1)
        mts=np.stack([c for i in range(s) for c in (ex[i][:,0]+ex[i][:,1],ex[i][:,2]+ex[i][:,3])],axis=1)
        sg=(1,-1)*s; w=0
    _,inv=np.unique(keys,axis=0,return_inverse=True); inv=inv.ravel()
    c=top.shape[1]
    order=np.lexsort(tuple(top[:,j] for j in range(c-1,-1,-1))+(inv,))  # group, then local tuple lexicographic
    bounds=np.flatnonzero(np.diff(inv[order]))+1
    R=[];C=[];V=[]; col=0
    for g in np.split(order,bounds):
        nloc,K=local_kernel(tuple(int(x) for x in mts[g[0]]),sg,w)
        assert nloc==len(g)
        if K.shape[0]==0: continue
        jj,ii=np.nonzero(K)
        R.append(g[ii]); C.append(col+jj); V.append(K[jj,ii]); col+=K.shape[0]
    if not R: return sps.csr_matrix((n,0),dtype=np.int64)
    return sps.csr_matrix((np.concatenate(V),(np.concatenate(R),np.concatenate(C))),shape=(n,col),dtype=np.int64)
# ---------- mod-P helpers ----------
def mulmod(A,B):
    # exact (A @ B) mod P for dense int64 arrays with entries in [0,P), via float64 split
    B1=(B>>10).astype(np.float64); B0=(B&1023).astype(np.float64); Af=A.astype(np.float64)
    hi=np.fmod(Af@B1,P).astype(np.int64); lo=np.fmod(Af@B0,P).astype(np.int64)
    return (hi*1024+lo)%P
def fill_nmod(Msp):
    Mc=Msp.tocoo(); A=flint.nmod_mat(Msp.shape[0],Msp.shape[1],P)
    for i,j,v in zip(Mc.row.tolist(),Mc.col.tolist(),(Mc.data%P).tolist()): A[i,j]=v
    return A
def nullspace_inplace(A):
    rows,cols=A.nrows(),A.ncols()
    _,rk=A.rref(inplace=True)
    piv=[]; j=0
    for i in range(rk):
        while int(A[i,j])==0: j+=1
        piv.append(j); j+=1
    ps=set(piv); free=[c for c in range(cols) if c not in ps]
    N=np.zeros((len(free),cols),dtype=np.int64)
    for t,f in enumerate(free):
        N[t,f]=1
        for i in range(rk):
            v=int(A[i,f])
            if v: N[t,piv[i]]=(P-v)%P
    return rk,N
def rank_dense(Tm):
    A=flint.nmod_mat(Tm.shape[0],Tm.shape[1],(Tm%P).ravel().tolist(),P); return A.rank()
# ---------- the test ----------
def run(kap,seed=5,compress=True):
    t0=time.time(); STATE['word']=kap; s=len(kap)
    m1,mU,mad=char_mults(kap,[(0,0),(1,1),(2,0)])
    W=Word(kap); n11=len(W.I((1,1)))
    log(f'{kap}: dimM={W.dimM}, dim M_(1,1)={n11}, M_(0,0)={len(W.I((0,0)))}, M_(2,0)={len(W.I((2,0)))}; '
        f'char (m1,mU,mad)=({m1},{mU},{mad}); 3mU={3*mU}, 5m1+mad={5*m1+mad}')
    if mU==0: log('  trivially injective'); return True,dict(m1=m1,mU=mU,mad=mad,rank=0)
    stage('root kernels'); ROOT=os.environ.get('ROOT','auto')
    cands=[1,2] if ROOT=='auto' else [int(ROOT)]
    Ks={r_:ker_root(W,r_) for r_ in cands}
    root=min(cands,key=lambda r_:Ks[r_].shape[1]); K2=Ks[root]; n2=K2.shape[1]; del Ks
    other,wo=((ea1,(2,0)) if root==2 else (ea2,(1,3)))
    stage('R K'); R1=W.op(other,(1,1),wo); A1=(R1@K2).tocsr(); A1.data%=P; A1.eliminate_zeros()
    A1=A1[np.flatnonzero(np.diff(A1.indptr))]
    rows=A1.shape[0]
    if compress and rows>n2+32:
        rng=np.random.default_rng(seed+1); d=n2+32; nnzc=4
        Sr=rng.integers(0,d,size=(rows*nnzc)); Sc=np.repeat(np.arange(rows),nnzc); Sv=rng.integers(1,P,size=rows*nnzc)
        S=sps.csr_matrix((Sv,(Sr,Sc)),shape=(d,rows),dtype=np.int64); S.data%=P
        A1=(S@A1).tocsr(); A1.data%=P
    log(f'  ker e_a{root}: {n2} of {n11} ({time.time()-t0:.1f}s); dense stage {A1.shape[0]}x{n2}, nnz={A1.nnz}')
    stage('fill'); t1=time.time(); A=fill_nmod(A1); del A1
    stage('rref'); rk,N=nullspace_inplace(A); del A
    log(f'  hw nullity {N.shape[0]} of {rows}x{n2} ({time.time()-t1:.1f}s); peakRSS={peak_gb():.2f}GB')
    if N.shape[0]!=mU:
        log(f'  nullity {N.shape[0]} != char m_U {mU}; retrying without row compression' if compress else '  MISMATCH')
        if compress: return run(kap,seed,compress=False)
        return False,dict(m1=m1,mU=mU,mad=mad,rank=-1)
    rng=np.random.default_rng(seed)
    a=rng.integers(1,P,size=(5,3,s)); b=rng.integers(1,P,size=(3,s))
    if os.environ.get('COEF')=='power':   # fixed slot-power coefficients (known to fail at (2,2,2))
        a=np.array([[[pow(i+1,r+1+5*c,P) for i in range(s)] for c in range(3)] for r in range(5)],dtype=np.int64)
        b=np.array([[pow(i+1,5+c,P) for i in range(s)] for c in range(3)],dtype=np.int64)
    n00=len(W.I((0,0))); n20=len(W.I((2,0)))
    dA=min(n00,m1+8); dB=min(n20,mad+8)
    stage('chain ops')
    Lops=[W.op(L,u,v) for L,(u,v) in zip((fa2,fa1,fa1,fa2),zip(Ywt[:-1],Ywt[1:]))]
    SA=[[W.op(Yd[k],Ywt[k],(0,0),slots=[i]) for k in range(5)] for i in range(s)]
    SB=[(W.op(Y[0],(1,-1),(2,0),slots=[i]),W.op(Y[1],(1,1),(2,0),slots=[i])) for i in range(s)]
    csz=int(max(16,min(mU,2e9/(8*5*max(n11,n00)))))
    for nnzc in (3,24):
        # sparse random projections of M_(0,0) -> F^dA and M_(2,0) -> F^dB (only lower rank; full rank still certifies)
        def proj(d,m):
            r_=rng.integers(0,d,size=m*nnzc); c_=np.repeat(np.arange(m),nnzc); v_=rng.integers(1,P,size=m*nnzc)
            M_=sps.csr_matrix((v_,(r_,c_)),shape=(d,m),dtype=np.int64); M_.data%=P; return M_
        PA=proj(dA,n00); PB=proj(dB,n20)
        QA=[[(PA@SA[i][k]).tocsr() for k in range(5)] for i in range(s)]
        QB=[((PB@SB[i][0]).tocsr(),(PB@SB[i][1]).tocsr()) for i in range(s)]
        for Q in [q for row in QA for q in row]+[q for pr in QB for q in pr]: Q.data%=P
        TA=np.zeros((3,5,dA,mU),dtype=np.int64); TB=np.zeros((3,dB,mU),dtype=np.int64)
        stage(f'slot maps (chunks of {csz}, proj nnz {nnzc})')
        for c0 in range(0,mU,csz):
            cs=slice(c0,min(mU,c0+csz))
            F=[(K2@N[cs].T)%P]
            for L in Lops: F.append((L@F[-1])%P)
            for i in range(s):
                x=QA[i][0]@F[0]
                for k in range(1,5): x=(x+QA[i][k]@F[k])%P
                y=(QB[i][0]@F[1]-QB[i][1]@F[0])%P
                for c in range(3):
                    for r in range(5): TA[c,r,:,cs]=(TA[c,r,:,cs]+int(a[r,c,i])*x)%P
                    TB[c,:,cs]=(TB[c,:,cs]+int(b[c,i])*y)%P
            del F
        Tm=np.concatenate([np.concatenate([TA[c,r] for r in range(5)]+[TB[c]],axis=0).T for c in range(3)],axis=0)
        del TA,TB
        stage('rank T'); rkT=rank_dense(Tm); del Tm
        if rkT==3*mU or nnzc==24: break
        log(f'  rank(T)={rkT} < {3*mU} with sparse projection nnz {nnzc}; retrying with denser projection')
    ok=rkT==3*mU
    log(f'  rank(T)={rkT} / 3mU={3*mU}: {"INJECTIVE" if ok else "NOT injective"}  (total {time.time()-t0:.0f}s, peakRSS={peak_gb():.2f}GB)')
    STATE['word']=None
    return ok,dict(m1=m1,mU=mU,mad=mad,rank=rkT)
if __name__=='__main__':
    threading.Thread(target=heartbeat,args=(float(os.environ.get('HB','60')),),daemon=True).start()
    for arg in sys.argv[1:]:
        run(tuple(int(x) for x in arg.split(',')))
