from collections import defaultdict
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from random import Random
import re
import time

BASE=Path('/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc')

def torus(B):
    A={(0,0):1}
    for z in B:
        n=abs(z); eps=1 if z>0 else -1; F=defaultdict(int)
        for r in range(n+1): F[r,r]+=1
        for r in range(n+1): F[r,n-r]+=eps
        O=defaultdict(int)
        for (i,j),x in A.items():
            for (u,v),y in F.items(): O[i+u,j+v]+=x*y
        A={k:v for k,v in O.items() if v}
    return A

def deformed(B):
    A={(0,0,0):1}
    for z in B:
        n=abs(z); eps=1 if z>0 else -1; F=defaultdict(int)
        for r in range(n+1): F[r,r,0]+=1
        for r in range(n//2+1):
            d=n-2*r; c=eps*(-1)**r*comb(n-r,r)
            for j in range(d+1): F[r+j,r+d-j,d]+=c*comb(d,j)
        O=defaultdict(int)
        for (i,j,e),x in A.items():
            for (u,v,f),y in F.items(): O[i+u,j+v,e+f]+=x*y
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

def CG(a,b):
    return range(abs(a-b),a+b+1,2)

class Core:
    def __init__(self,B):
        self.T=sum(map(abs,B))
        self.M=torus(B)
        self.hcache={}
        self.jcache={}
    def m(self,i,j):
        if min(i,j)<0 or max(i,j)>self.T:
            return 0
        return self.M.get((i,j),0)
    def H(self,i,j):
        if min(i,j)<0:
            return 0
        if i<j:
            i,j=j,i
        key=(i,j)
        if key not in self.hcache:
            self.hcache[key]=self.m(i,j)-self.m(i+1,j-1)
        return self.hcache[key]
    def P(self,j):
        return self.H(j,j)
    def J(self,k,i,j):
        if min(i,j)<0:
            return Q(0)
        if i<j:
            i,j=j,i
        key=(k,i,j)
        if key not in self.jcache:
            alpha=Q(k-i+j,2)
            beta=Q(k+i-j,2)
            omega=1/(beta+1)
            value=Q(0)
            for t in range(j+1):
                value+=omega*self.H(i+t,j-t)
                omega*=(alpha-t)/(beta+t+2)
            self.jcache[key]=value
        return self.jcache[key]
    def blocks(self,B,p,R):
        W=sum(map(abs,B))
        delta=(W-p)//2
        za,zb=sorted(R,key=abs,reverse=True)
        a,b=abs(za),abs(zb)
        ea=1 if za>0 else -1
        eb=1 if zb>0 else -1
        r=delta-a
        s=delta-b
        ell=delta-a-b-1
        h=delta-(a+b)//2
        S=sum(self.P(t) for t in range(s,delta+1))-sum(self.P(t) for t in range(ell,r))
        Z=ea*eb*(self.m(s,r)-self.m(delta+1,ell)-self.m(s-1,r-1)+self.m(delta,ell-1))
        child=self.P(h)-self.P(h-1)
        A=eb*self.H(delta,s)+ea*self.H(delta,r)
        B0=eb*self.H(r-1,ell)+ea*self.H(s-1,ell)
        Q0=S+Z-child-B0
        return delta,a,b,ea,eb,r,s,ell,S,Z,child,A,B0,Q0
    def energy(self,delta,a,b,ea,eb,r,s,k):
        V=(b+1)**2*self.J(2*b-k,s,s)+(a+1)**2*self.J(2*a-k,r,r)
        V+=2*ea*eb*(a+1)*(b+1)*self.J(a+b-k,r,s)
        return self.J(k,delta,delta)*V

def top_pair(B):
    pairs=[]
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if (B[i]-B[j])%2==0:
                pairs.append((abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])),B[i],B[j]))
    if not pairs:
        raise ValueError('no same-parity pair')
    q=max(pairs,key=lambda z:(z[0],z[1]))
    return q[2],q[3]

def parse_file(path,kind):
    out=[]
    for line in path.read_text().splitlines():
        if not line.startswith(kind+' '):
            continue
        B=tuple(map(int,re.search(r'\bB=(.*?)(?:  |$)',line)[1].split()))
        p=abs(int(re.search(r'\bp=(-?\d+)',line)[1]))
        W=sum(map(abs,B))
        delta=(W-p)//2
        assert (W-p)%2==0 and p>=max(6,max(map(abs,B)))
        assert delta>=8 and max(map(abs,B))<=delta and all(-z not in B for z in B)
        out.append((p,B))
    return out

def moment_screen():
    rng=Random(161103)
    profiles=[(),(1,),(-2,),(-3,2),(-3,-3,1,1),(-1,2,3,-4)]
    profiles += [tuple(rng.choice((-1,1))*rng.randrange(1,5)
                       for _ in range(rng.randrange(1,5))) for _ in range(30)]
    checks=0
    for B in profiles:
        W=sum(map(abs,B))
        M=torus(B)
        V=deformed(B)
        def h(i,j):
            if min(i,j)<0:
                return 0
            if i<j:
                i,j=j,i
            return M.get((i,j),0)-M.get((i+1,j-1),0)
        for k in range(7):
            direct=defaultdict(Q)
            for (i,j,e),v in V.items():
                direct[i,j]+=Q(2*v,k+e+2)
            for i in range(W+1):
                for j in range(i+1):
                    alpha=Q(k-i+j,2)
                    beta=Q(k+i-j,2)
                    omega=1/(beta+1)
                    value=Q(0)
                    for t in range(j+1):
                        value+=omega*h(i+t,j-t)
                        omega*=(alpha-t)/(beta+t+2)
                    assert value==direct[i,j],(B,k,i,j)
                    checks+=1
    return len(profiles),checks

def random_screen():
    rng=Random(161103)
    checked=passed=cschecks=0
    for _ in range(900):
        vals=rng.sample(range(1,11),rng.randrange(4,9))
        B=tuple(v if rng.randrange(2) else -v for v in vals)
        W=sum(map(abs,B))
        mx=max(map(abs,B))
        low=max(8,mx)
        high=(W-max(6,mx))//2
        if high<low:
            continue
        delta=rng.randrange(low,high+1)
        p=W-2*delta
        pairs=[(B[i],B[j]) for i in range(len(B)) for j in range(i+1,len(B))
               if (B[i]-B[j])%2==0]
        if not pairs:
            continue
        R=rng.choice(pairs)
        C=list(B)
        C.remove(R[0])
        C.remove(R[1])
        C=tuple(C)
        core=Core(C)
        vals=core.blocks(B,p,R)
        d,a,b,ea,eb,r,s,ell,S,Z,child,A,B0,Q0=vals
        G=char_table(C)
        parent=char_table(B).get((p,0),0)
        child_direct=G.get((p,0),0)
        X=eb*sum(G.get((u,b),0) for u in CG(p,a))
        Y=ea*sum(G.get((u,a),0) for u in CG(p,b))
        Sdirect=sum(G.get((u,0),0) for t in CG(a,b) for u in CG(p,t))
        Zdirect=ea*eb*sum(G.get((p,t),0) for t in CG(a,b))
        assert S==Sdirect and Z==Zdirect and child==child_direct
        assert X+Y==A-B0 and parent-child==Q0+A
        energies=[core.energy(d,a,b,ea,eb,r,s,k) for k in range(2*b+1)]
        assert all(e>=0 for e in energies)
        assert all(A*A<=e for e in energies)
        cschecks+=len(energies)
        checked+=1
        if Q0>=0 and Q0*Q0>=min(energies):
            passed+=1
            assert parent>=child,(B,p,R,parent,child)
    return checked,passed,cschecks

def census(rows,label,all_k):
    groups=defaultdict(list)
    for p,B in rows:
        R=top_pair(B)
        C=list(B)
        C.remove(R[0])
        C.remove(R[1])
        groups[tuple(sorted(C))].append((p,B,R))
    certs=rejects=0
    khist=defaultdict(int)
    failed=[]
    start=time.time()
    for n,(C,items) in enumerate(groups.items(),1):
        core=Core(C)
        for p,B,R in items:
            data=core.blocks(B,p,R)
            d,a,b,ea,eb,r,s,ell,S,Z,child,A,B0,Q0=data
            if all_k:
                ks=range(2*b+1)
            else:
                ks=sorted({0,b-b%2,b+b%2})
            energies=[(k,core.energy(d,a,b,ea,eb,r,s,k)) for k in ks]
            assert all(E>=0 for k,E in energies)
            bestk,best=min(energies,key=lambda z:z[1])
            if Q0>=0 and Q0*Q0>=best:
                certs+=1
                khist[bestk]+=1
            else:
                rejects+=1
                failed.append((B,p,R,Q0,bestk,best))
        if n%1000==0 and len(groups)>2000:
            print(f'{label}: child tables {n}/{len(groups)}, rows {certs+rejects}/{len(rows)}',flush=True)
    print(label,'rows=',len(rows),'children=',len(groups),'certified=',certs,'rejected=',rejects,
          'best-k=',dict(khist),'seconds=',round(time.time()-start,1),flush=True)
    return certs,rejects,failed

profiles,identity_count=moment_screen()
print('direct moment profiles/identities=',profiles,identity_count)
print('random blocks, CB descents, C-S pairs=',random_screen())

no_flip=parse_file(BASE/'fx3b_48.log','NOFLIP')
controls=[]
for name in ('fx5_40.log','tpg_one.log','tpg_all.log'):
    controls+=parse_file(BASE/name,'TPFAIL')
assert len(no_flip)==33487 and len(controls)==190

nf_result=census(no_flip,'no-flip candidate k',False)
control_result=census(controls,'TopPair controls, all k',True)
assert nf_result[:2]==(33487,0) and control_result[:2]==(0,190)
print('ALL EXACT CHECKS PASS')
