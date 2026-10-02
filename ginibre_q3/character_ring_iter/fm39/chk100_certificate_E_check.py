from collections import defaultdict
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from random import Random
import re,time

LOG=Path('/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc')

def torus(B):
    A={(0,0):1}
    for z in B:
        n=abs(z); eps=1 if z>0 else -1; Fct=defaultdict(int)
        for r in range(n+1): Fct[r,r]+=1
        for r in range(n+1): Fct[r,n-r]+=eps
        O=defaultdict(int)
        for (i,j),x in A.items():
            for (u,v),y in Fct.items(): O[i+u,j+v]+=x*y
        A={k:v for k,v in O.items() if v}
    return A

def cpoly(B):
    A={(0,0,0):1}
    for z in B:
        n=abs(z); eps=1 if z>0 else -1; Fct=defaultdict(int)
        for r in range(n+1): Fct[r,r,0]+=1
        for r in range(n//2+1):
            d=n-2*r; c=eps*(-1)**r*comb(n-r,r)
            for j in range(d+1): Fct[r+j,r+d-j,d]+=c*comb(d,j)
        O=defaultdict(int)
        for (i,j,e),x in A.items():
            for (u,v,f),y in Fct.items(): O[i+u,j+v,e+f]+=x*y
        A={k:v for k,v in O.items() if v}
    return A

def char_table(B):
    A={(0,0):1}
    for z in B:
        n=abs(z); eps=1 if z>0 else -1; O=defaultdict(int)
        for (i,j),v in A.items():
            for k in range(abs(i-n),i+n+1,2): O[k,j]+=v
            for k in range(abs(j-n),j+n+1,2): O[i,k]+=eps*v
        A={k:v for k,v in O.items() if v}
    return A

class Core:
    def __init__(self,B):
        self.T=sum(map(abs,B)); self.M=torus(B); self.hc={}; self.jc={}
    def m(self,i,j):
        if min(i,j)<0 or max(i,j)>self.T: return 0
        return self.M.get((i,j),0)
    def h(self,i,j):
        if min(i,j)<0: return 0
        if i<j: i,j=j,i
        key=(i,j)
        if key not in self.hc: self.hc[key]=self.m(i,j)-self.m(i+1,j-1)
        return self.hc[key]
    def P(self,j): return self.h(j,j)
    def J(self,k,i,j):
        if min(i,j)<0: return Q(0)
        if i<j: i,j=j,i
        key=(k,i,j)
        if key in self.jc: return self.jc[key]
        t=i-j; beta=Q(k+t,2); alpha=Q(k-t,2)
        omega=1/(beta+1); value=Q(0)
        for s in range(j+1):
            value+=omega*self.h(i+s,j-s)
            omega*= (alpha-s)/(beta+s+2)
        self.jc[key]=value
        return value
    def certificate(self,B,p,R):
        W=sum(map(abs,B)); delta=(W-p)//2
        za,zb=sorted(R,key=abs,reverse=True); a,b=abs(za),abs(zb)
        ea=1 if za>0 else -1; eb=1 if zb>0 else -1
        r=delta-a; s=delta-b; ell=delta-a-b-1; h=delta-(a+b)//2
        S=sum(self.P(j) for j in range(s,delta+1))-sum(self.P(j) for j in range(ell,r))
        X=eb*(self.h(delta,s)-self.h(r-1,ell))
        Y=ea*(self.h(delta,r)-self.h(s-1,ell))
        Z=ea*eb*(self.m(s,r)-self.m(delta+1,ell)-self.m(s-1,r-1)+self.m(delta,ell-1))
        gc=self.P(h)-self.P(h-1)
        U=(b+1)**2*self.J(2*b,s,s)+(a+1)**2*self.J(2*a,r,r)+2*ea*eb*(a+1)*(b+1)*self.J(a+b,r,s)
        V=Q(0)
        if ell>=0:
            V=(b+1)**2*self.J(2*b,r-1,r-1)+(a+1)**2*self.J(2*a,s-1,s-1)+2*ea*eb*(a+1)*(b+1)*self.J(a+b,r-1,s-1)
        q=S+Z-gc; ep=self.P(delta)*U; em=self.P(ell)*V
        D=Q(q*q)-ep-em
        E=(q>=0 and D>=0 and D*D>=4*ep*em)
        return {'E':E,'q':q,'D':D,'ep':ep,'em':em,'S':S,'X':X,'Y':Y,'Z':Z,'gc':gc,'blocks':S+X+Y+Z}

def top_pair(B):
    choices=[]
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if (B[i]-B[j])%2==0:
                choices.append((abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])),B[i],B[j]))
    if not choices: raise ValueError('no same-parity pair')
    q=max(choices,key=lambda t:(t[0],t[1]))
    return q[2],q[3]

def parse(path,kind):
    rows=[]
    for line in path.read_text().splitlines():
        if not line.startswith(kind+' '): continue
        B=tuple(map(int,re.search(r'\bB=(.*?)(?:  |$)',line)[1].split()))
        p=abs(int(re.search(r'\bp=(-?\d+)',line)[1]))
        W=sum(map(abs,B)); delta=(W-p)//2
        assert (W-p)%2==0 and p>=max(6,max(map(abs,B)))
        assert delta>=8 and max(map(abs,B))<=delta and all(-z not in B for z in B)
        logged=re.search(r'TopPair\((-?\d+),(-?\d+)\)',line)
        if logged:
            chosen=top_pair(B)
            assert sorted(chosen)==sorted((int(logged[1]),int(logged[2])))
        rows.append((p,B))
    return rows

def census(rows,label):
    groups=defaultdict(list)
    for p,B in rows:
        R=top_pair(B); C=list(B); C.remove(R[0]); C.remove(R[1])
        groups[tuple(sorted(C))].append((p,B,R))
    good=bad=0; bad_rows=[]; started=time.time()
    for ix,(C,records) in enumerate(groups.items(),1):
        core=Core(C)
        for p,B,R in records:
            z=core.certificate(B,p,R)
            if z['E']: good+=1
            else: bad+=1; bad_rows.append((B,p,R,z))
        if ix%1000==0 and len(groups)>2000:
            print(f'{label}: child tables {ix}/{len(groups)}, rows {good+bad}/{len(rows)}',flush=True)
    print(f'{label}: rows={len(rows)} child-tables={len(groups)} E_pass={good} E_fail={bad}',flush=True)
    return good,bad,bad_rows

def factor_screen():
    cs=(Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)); count=0
    for n in range(13):
        for eps in (-1,1):
            for c in cs:
                mat=[[Q(0) for _ in range(n+1)] for _ in range(n+1)]
                for r in range(n+1): mat[r][r]+=1
                for r in range(n//2+1):
                    d=n-2*r; coef=eps*(-1)**r*comb(n-r,r)*c**d
                    for j in range(d+1): mat[r+j][r+d-j]+=coef*comb(d,j)
                ds=[]
                for j in range(n+1):
                    dj=sum((Q((-1)**l*comb(j,l)*comb(n-j,l))*c**(n-2*l)*(1-c*c)**l
                            for l in range(min(j,n-j)+1)),Q(0))
                    ds.append(dj); assert abs(dj)<=1
                for i in range(n+1):
                    for j in range(n+1):
                        assert mat[i][j]==(Q(i==j)+eps*(ds[i] if i+j==n else 0))
                count+=1
    return count

def exact_psd(A):
    A=[list(map(Q,row)) for row in A]; n=len(A)
    while n:
        if any(A[i][i]<0 for i in range(n)): return False
        pivot=next((i for i in range(n) if A[i][i]>0),None)
        if pivot is None: return all(A[i][j]==0 for i in range(n) for j in range(n))
        if pivot:
            A[pivot],A[0]=A[0],A[pivot]
            for row in A: row[pivot],row[0]=row[0],row[pivot]
        p=A[0][0]
        A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,n)] for i in range(1,n)]
        n-=1
    return True

def kernel_psd_screen():
    rng=Random(15699); cs=(Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)); checks=0
    for _ in range(30):
        B=tuple(rng.choice((-1,1))*rng.randrange(1,4) for _ in range(rng.randrange(0,4)))
        T=sum(map(abs,B)); A=cpoly(B)
        for c in cs:
            M=[[sum((Q(v)*c**e for (i,j,e),v in A.items() if (i,j)==(u,vv)),Q(0))
                for vv in range(T+1)] for u in range(T+1)]
            assert all(M[i][j]==M[j][i] for i in range(T+1) for j in range(T+1))
            assert exact_psd(M),(B,c)
            checks+=1
    return checks

def moment_screen():
    rng=Random(156100); profiles=[(),(1,),(-2,),(-3,2),(-3,-3,1,1),(-1,2,3,-4)]
    profiles += [tuple(rng.choice((-1,1))*rng.randrange(1,5)
                       for _ in range(rng.randrange(1,5)) for __ in [0])
                 for _ in range(30)]
    checks=0
    for B in profiles:
        T=sum(map(abs,B)); A=cpoly(B); M=torus(B)
        def m(i,j): return M.get((i,j),0) if min(i,j)>=0 and max(i,j)<=T else 0
        def H(i,j):
            if min(i,j)<0: return 0
            if i<j: i,j=j,i
            return m(i,j)-m(i+1,j-1)
        for i in range(T+1):
            for j in range(T+1):
                assert sum(v for (x,y,e),v in A.items() if (x,y)==(i,j))==m(i,j)
        def formula(k,i,j):
            if min(i,j)<0: return Q(0)
            if i<j: i,j=j,i
            t=i-j; beta=Q(k+t,2); alpha=Q(k-t,2)
            omega=1/(beta+1); value=Q(0)
            for s in range(j+1):
                value+=omega*H(i+s,j-s); omega*=(alpha-s)/(beta+s+2)
            return value
        for k in range(5):
            direct=defaultdict(Q)
            for (i,j,e),v in A.items(): direct[i,j]+=Q(2*v,k+e+2)
            for i in range(T+1):
                for j in range(T+1):
                    assert direct[i,j]==formula(k,i,j),(B,k,i,j)
                    checks+=1
    return len(profiles),checks

def random_descent_screen():
    rng=Random(156200); checked=passed=0
    for _ in range(800):
        vals=rng.sample(range(1,10),rng.randrange(4,9))
        B=tuple(v if rng.randrange(2) else -v for v in vals)
        W=sum(map(abs,B)); mx=max(map(abs,B)); lo=max(8,mx); hi=(W-max(6,mx))//2
        if hi<lo: continue
        delta=rng.randrange(lo,hi+1); p=W-2*delta
        pairs=[(B[i],B[j]) for i in range(len(B)) for j in range(i+1,len(B))
               if (B[i]-B[j])%2==0]
        if not pairs: continue
        R=rng.choice(pairs); C=list(B); C.remove(R[0]); C.remove(R[1])
        z=Core(C).certificate(B,p,R)
        parent=char_table(B).get((p,0),0); child=char_table(C).get((p,0),0)
        assert z['blocks']==parent and z['gc']==child
        checked+=1
        if z['E']:
            passed+=1
            assert parent>=child,(B,p,R,parent,child,z)
    return checked,passed

# Direct exact integration of the c-polynomial, independent of the J recurrence.
def direct_E(B,p,R):
    W=sum(map(abs,B)); delta=(W-p)//2; C=list(B)
    C.remove(R[0]); C.remove(R[1]); A=cpoly(C)
    def m(i,j): return sum(v for (x,y,e),v in A.items() if (x,y)==(i,j))
    def H(i,j):
        if min(i,j)<0: return 0
        if i<j: i,j=j,i
        return m(i,j)-m(i+1,j-1)
    def P(i): return H(i,i)
    def J(k,i,j):
        if min(i,j)<0: return Q(0)
        return sum((Q(2*v,k+e+2) for (x,y,e),v in A.items() if (x,y)==(i,j)),Q(0))
    za,zb=sorted(R,key=abs,reverse=True); a,b=abs(za),abs(zb)
    ea=1 if za>0 else -1; eb=1 if zb>0 else -1
    r=delta-a; s=delta-b; ell=delta-a-b-1; h=delta-(a+b)//2
    S=sum(P(j) for j in range(s,delta+1))-sum(P(j) for j in range(ell,r))
    Z=ea*eb*(m(s,r)-m(delta+1,ell)-m(s-1,r-1)+m(delta,ell-1))
    gc=char_table(C).get((p,0),0)
    U=(b+1)**2*J(2*b,s,s)+(a+1)**2*J(2*a,r,r)+2*ea*eb*(a+1)*(b+1)*J(a+b,r,s)
    V=Q(0)
    if ell>=0:
        V=(b+1)**2*J(2*b,r-1,r-1)+(a+1)**2*J(2*a,s-1,s-1)+2*ea*eb*(a+1)*(b+1)*J(a+b,r-1,s-1)
    q=S+Z-gc; ep=P(delta)*U; em=P(ell)*V; D=Q(q*q)-ep-em
    E=(q>=0 and D>=0 and D*D>=4*ep*em)
    return q,D,ep,em,char_table(B).get((p,0),0),gc,E

print('factor decomposition cases=',factor_screen())
print('whole-kernel exact PSD samples=',kernel_psd_screen())
print('moment profiles, exact identities=',moment_screen())
print('random residual pair checks (E-passing)=',random_descent_screen())

rows_nf40=parse(LOG/'fx5_40.log','NOFLIP')
rows_tf40=parse(LOG/'fx5_40.log','TPFAIL')
rows_nf48=parse(LOG/'fx3b_48.log','NOFLIP')
assert len(rows_nf40)==5430 and len(rows_tf40)==126 and len(rows_nf48)==33487
r1=census(rows_nf40,'W<=40 no-flip')
r2=census(rows_tf40,'W<=40 TopPair failures')
r3=census(rows_nf48,'W<=48 no-flip')
assert r1[:2]==(5430,0) and r2[:2]==(0,126) and r3[:2]==(33485,2)
print('raw W<=48 failing records:')
for B,p,R,z in r3[2]:
    direct=direct_E(B,p,R)
    assert direct[:4]==(z['q'],z['D'],z['ep'],z['em'])
    assert direct[4:6]==(4075371,160741) and not direct[6]
    print('B=',B,'p=',p,'R=',R,'q,D,E+,E-=',direct[:4],'parent,child=',direct[4:6])

extra=[
    (6,(1,1,-2,-2,3,3,3,3,-4,-4,-4,5,5,5,5)),
    (6,(-1,-1,-2,-2,-3,-3,-3,-3,-4,-4,-4,-5,-5,-5,-5)),
    (8,(1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)),
    (8,(-1,-1,-2,-3,-3,-3,-3,-4,-4,-4,-5,-5,-5,-5,8)),
]
for p,B in extra:
    R=top_pair(B); z=direct_E(B,p,R)
    assert sum(map(abs,B))>48 and not z[6]
    print('outside-range E failure: W=',sum(map(abs,B)),'p=',p,'B=',B,'R=',R,'q,D,E+,E-=',z[:4])
print('ALL EXACT CHECKS PASS')
