# Mixed q_ij Fock space, 2 colours: q11=q22=q, q12=q21=s. Exact rationals via Fraction.
# F(T) = < prod_i (H_{n_i}(X_1|q) + eps_i H_{n_i}(X_2|q)) Omega, Omega >, operators applied right-to-left.
import itertools, random, sys
from fractions import Fraction as Fr
def make(q,s):
    Q={(0,0):q,(1,1):q,(0,1):s,(1,0):s}
    def ann(i,vec):  # vec: dict word(tuple)->coef
        out={}
        for w,c in vec.items():
            coef=c
            for k,col in enumerate(w):
                if col==i:
                    nw=w[:k]+w[k+1:]
                    out[nw]=out.get(nw,0)+coef
                coef=coef*Q[(i,col)]
                if coef==0: break
        return out
    def cre(i,vec):
        return {(i,)+w:c for w,c in vec.items()}
    def X(i,vec,maxdeg):
        a=ann(i,vec); b={w:c for w,c in cre(i,vec).items() if len(w)<=maxdeg}
        for w,c in b.items(): a[w]=a.get(w,0)+c
        return {w:c for w,c in a.items() if c!=0}
    def qint(n): return sum(q**k for k in range(n))
    def H(n,i,vec,maxdeg):
        # H_{k+1} = X H_k - [k] H_{k-1}, applied to vec
        prev={}; cur=dict(vec)
        for k in range(n):
            nxt=X(i,cur,maxdeg)
            if k>=1:
                for w,c in prev.items(): nxt[w]=nxt.get(w,0)-qint(k)*c
            prev,cur=cur,{w:c for w,c in nxt.items() if c!=0}
        return cur
    return H
def F(labels,T,q,s):
    H=make(q,s); L=len(labels)
    vec={():Fr(1)}
    suffix=[0]*(L+1)
    for i in range(L-1,-1,-1): suffix[i]=suffix[i+1]+labels[i]
    for i in range(L-1,-1,-1):
        maxdeg=min(sum(labels[:i]), suffix[i]) # degree after applying block i must be <= what remaining blocks (0..i-1) can remove
        maxdeg=sum(labels[:i])
        e=-1 if T>>i&1 else 1
        a=H(labels[i],0,vec,maxdeg+labels[i]); b=H(labels[i],1,vec,maxdeg+labels[i])
        new={}
        for w,c in a.items():
            if len(w)<=maxdeg: new[w]=new.get(w,0)+c
        for w,c in b.items():
            if len(w)<=maxdeg: new[w]=new.get(w,0)+e*c
        vec={w:c for w,c in new.items() if c!=0}
    return vec.get((),0)
if __name__=="__main__":
    # self-check against the enumeration in qs_family_screen.py
    src=open("qs_family_screen.py").read().split("import random")[0]
    ns={}; exec(src, ns)
    from fractions import Fraction as Fr
    bad=0
    for lab in [(1,1,1,1),(1,5,2,2),(1,1,2,2),(1,2,3)]:
        for T in range(1<<len(lab)):
            if bin(T).count("1")%2: continue
            p=ns["Fqs"](lab,T)
            for q,s_ in [(Fr(0),Fr(1)),(Fr(1,3),Fr(1,2)),(Fr(1,2),Fr(-1,5))]:
                bad+= F(lab,T,q,s_)!=ns["ev"](p,q,s_)
    print("self-check mismatches:",bad)
