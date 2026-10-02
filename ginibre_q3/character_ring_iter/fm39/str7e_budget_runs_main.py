# Interior height budget (FM-STR7e) for one list and pair; exact Python ints.
import sys, time
from collections import defaultdict
def cg(a,b): return range(abs(a-b),a+b+1,2)
def mul_factor(d,z):
    n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
    for (a,b),v in d.items():
        for c in cg(a,n): q[c,b]+=v
        for c in cg(b,n): q[a,c]+=e*v
    return {k:v for k,v in q.items() if v}
def table(word):
    d={(0,0):1}
    for z in sorted(word,key=abs,reverse=True): d=mul_factor(d,z)
    return d
def mul_xx(d,a,b,coef_y):  # multiply by U_a(x)U_b(x) + coef_y U_a(y)U_b(y)
    q=defaultdict(int)
    for (r,s),v in d.items():
        for c in cg(a,b):
            for r2 in cg(r,c): q[r2,s]+=v
            if coef_y:
                for s2 in cg(s,c): q[r,s2]+=coef_y*v
    return q
def mul_xy(d,a,b,ea,eb):  # multiply by eb U_a(x)U_b(y) + ea U_a(y)U_b(x)
    q=defaultdict(int)
    for (r,s),v in d.items():
        for r2 in cg(r,a):
            for s2 in cg(s,b): q[r2,s2]+=eb*v
        for r2 in cg(r,b):
            for s2 in cg(s,a): q[r2,s2]+=ea*v
    return q
def interior_cut(C):
    A=[];B=[];wa=wb=0
    for z in sorted(C,key=lambda z:(-abs(z),z)):
        if wa<=wb: A.append(z); wa+=abs(z)
        else: B.append(z); wb+=abs(z)
    if wa>wb: A,B=B,A
    return A,B
def budget(L,iu,iv):
    zu,zv=L[iu],L[iv]; a,b=abs(zu),abs(zv); eu=1 if zu>0 else -1; ev=1 if zv>0 else -1
    C=[z for i,z in enumerate(L) if i not in (iu,iv)]
    A,B=interior_cut(C)
    rA=table(A); FB=table(B)
    Gp=mul_xx(FB,a,b,eu*ev); Gm=mul_xy(FB,a,b,eu,ev)
    P=defaultdict(int); M=defaultdict(int)
    for (r,s),v in rA.items():
        P[r+s]+=v*Gp.get((r,s),0); M[r+s]+=v*Gm.get((r,s),0)
    run=0; mn=None; phi=0
    for t in sorted(set(P)|set(M)):
        run+=P[t]+min(M[t],0); phi+=P[t]+M[t]
        mn=run if mn is None or run<mn else mn
    return mn,phi,run
def toppair(B):
    best=None
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if (abs(B[i])-abs(B[j]))%2: continue
            key=(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])))
            if best is None or key>best[0]: best=(key,i,j)
    return best[1],best[2]
if __name__=="__main__":
    for k in range(int(sys.argv[1]),int(sys.argv[2])+1):
        sig=(-1)**k
        for p in range(k+1,k+4):
            if (k*(k+1)//2+p)%2: continue
            B=[-n for n in range(1,k+1)]; L=B+[sig*p]
            i,j=toppair(B); t0=time.time()
            mn,phi,fin=budget(L,i,j)
            print(time.strftime('%H:%M:%S'),f"k={k} p={sig*p} TopPair=({L[i]},{L[j]}) Phi={phi} minB/Phi={mn/phi:.4f} Binf/Phi={fin/phi:.4f} {'OK' if mn>=0 else 'FAIL'} {time.time()-t0:.1f}s",flush=True)
