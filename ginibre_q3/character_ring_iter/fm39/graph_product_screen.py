# Graph-product EVEN form: colours x_0..x_{N-1}, y_0..y_{N-1}; q=0 within/among families and on
# non-edges; q=1 on edges (x_a,y_b) in G.  x = N^-1/2 sum x_a, y = N^-1/2 sum y_b (exactly semicircular).
# F_G(T) = < prod_i (U_{n_i}(x) + eps_i U_{n_i}(y)) Omega, Omega >.  Exact arithmetic with sqrt(N) avoided
# by using unnormalized sums and rescaling: X = sum_a X(e_a) has variance N; use U_n(X/sqrt N) via
# the recursion U_{k+1}(t) = t U_k - U_{k-1} applied as  sqrtN*U_{k+1} = X U_k - sqrtN U_{k-1}... keep exact
# by working with V_k = N^{k/2} U_k(X/sqrtN): V_{k+1} = X V_k - N V_{k-1}; then prod over blocks rescales by N^{-D/2}.
import itertools, sys
from fractions import Fraction as Fr
def evaluator(N, edges):
    C=2*N
    def qq(i,j):
        if i<N and j>=N: return 1 if (i,j-N) in edges else 0
        if j<N and i>=N: return 1 if (j,i-N) in edges else 0
        return 0
    Q=[[qq(i,j) for j in range(C)] for i in range(C)]
    def ann(i,vec):
        out={}
        for w,c in vec.items():
            coef=c
            for k,col in enumerate(w):
                if col==i:
                    nw=w[:k]+w[k+1:]; out[nw]=out.get(nw,0)+coef
                coef=coef*Q[i][col]
                if coef==0: break
        return out
    def Xfam(fam,vec,maxdeg):
        cols=range(N) if fam==0 else range(N,2*N)
        out={}
        for i in cols:
            for w,c in ann(i,vec).items(): out[w]=out.get(w,0)+c
            for w,c in vec.items():
                if len(w)+1<=maxdeg:
                    nw=(i,)+w; out[nw]=out.get(nw,0)+c
        return {w:c for w,c in out.items() if c}
    def V(n,fam,vec,maxdeg):
        prev={}; cur=dict(vec)
        for k in range(n):
            nxt=Xfam(fam,cur,maxdeg)
            if k>=1:
                for w,c in prev.items(): nxt[w]=nxt.get(w,0)-N*c
            prev,cur=cur,{w:c for w,c in nxt.items() if c}
        return cur
    def F(labels,T):
        vec={():Fr(1)}; L=len(labels)
        for i in range(L-1,-1,-1):
            maxdeg=sum(labels[:i]); e=-1 if T>>i&1 else 1
            a=V(labels[i],0,vec,maxdeg+labels[i]); b=V(labels[i],1,vec,maxdeg+labels[i])
            new={}
            for w,c in a.items():
                if len(w)<=maxdeg: new[w]=new.get(w,0)+c
            for w,c in b.items():
                if len(w)<=maxdeg: new[w]=new.get(w,0)+e*c
            vec={w:c for w,c in new.items() if c}
        return vec.get((),0)/Fr(N)**(sum(labels)//2)
    return F
if __name__=="__main__":
    N=int(sys.argv[1]); maxD=int(sys.argv[2])
    allpairs=[(a,b) for a in range(N) for b in range(N)]
    # graphs up to relabelling: just enumerate all subsets (small N)
    graphs=[]
    for mask in range(1<<len(allpairs)):
        graphs.append(frozenset(p for k,p in enumerate(allpairs) if mask>>k&1))
    lists=[lab for L in range(2,7) for lab in itertools.combinations_with_replacement(range(1,4),L) if sum(lab)%2==0 and sum(lab)<=maxD]
    # sanity: complete graph = FM3 (classical), empty = free
    from math import comb
    neg=0; tot=0; worst=(1e9,None)
    seen=set()
    for G in graphs:
        # canonical form under S_N x S_N to skip isomorphic graphs
        can=min(tuple(sorted((pa[a],pb[b]) for a,b in G)) for pa in itertools.permutations(range(N)) for pb in itertools.permutations(range(N)))
        if can in seen: continue
        seen.add(can)
        F=evaluator(N,G)
        for lab in lists:
            for T in range(1<<len(lab)):
                if bin(T).count('1')%2: continue
                v=F(lab,T); tot+=1
                if v<worst[0]: worst=(v,(sorted(G),lab,T))
                if v<0: neg+=1
        print("graph",sorted(G),"done; profiles so far",tot,"neg",neg,flush=True)
    print("N",N,"graphs",len(seen),"profiles",tot,"negative",neg,"worst",worst)
