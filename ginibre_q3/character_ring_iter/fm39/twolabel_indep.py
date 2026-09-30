# Independent check of the FM-MECH7 two-label reduction with the main agent's phi code, and of (E) on/off the strip.
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
from math import comb
def cvec(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return c
def TW(a,e,p,q):
    if (a+e+p+q)%2: return 0,0
    p,q=max(p,q),min(p,q); c=cvec(a,e); N=a+e
    C=lambda k: c[k] if 0<=k<=N else 0
    B=lambda k: C(k-1)+C(k+1); D=lambda k: C(k)**2-C(k-1)*C(k+1)
    j=(N+p-q)//2; i=(N+p+q)//2+1
    return D(j)-D(i), B(i)*C(j)-C(i)*B(j)
bad=0; n=0
for r in range(1,5):
    for a in range(0,9):
        for s_ in range(1,6):
            for t_ in range(1,6):
                # HH: e=2r-2, labels s+1,t+1
                T,W=TW(a,2*r-2,s_+1,t_+1); v=phi(r,mul(h(s_),h(t_)),a); n+=1
                if v!=T-W: bad+=1
                T,W=TW(a,2*r,s_,t_); v=phi(r,mul(shat(s_),shat(t_)),a) if s_>=1 else None; n+=1
                if v!=T+W: bad+=1
                p,q=s_+1,t_; T,W=TW(a,2*r-1,p,q); v=phi(r,mul(h(s_),shat(t_)),a); n+=1
                if v!=(T+W if p>=q else T-W): bad+=1
print('reduction vs main-agent phi:',n,'cases, mismatches',bad)
# (E) on the strip |a-e|<=1 exhaustively e<=40, and off-strip evidence a,e<=24
def Echeck(a,e):
    c=cvec(a,e); N=a+e
    C=lambda k: c[k] if 0<=k<=N else 0
    B=lambda k: C(k-1)+C(k+1); D=lambda k: C(k)**2-C(k-1)*C(k+1)
    f=0; m=0
    for j in range((N+1)//2,N+2):
        for i in range(j+1,N+3):
            m+=1
            if D(j)-D(i) < abs(B(i)*C(j)-C(i)*B(j)): f+=1
    return f,m
fs=ms=0
for e in range(0,41):
    for a in (e-1,e,e+1):
        if a<0: continue
        f,m=Echeck(a,e); fs+=f; ms+=m
print('strip |a-e|<=1, e<=40:',ms,'pairs, failures',fs)
fs=ms=0
for a in range(0,25):
    for e in range(0,25):
        f,m=Echeck(a,e); fs+=f; ms+=m
print('all a,e<=24 (evidence for open region):',ms,'pairs, failures',fs)
