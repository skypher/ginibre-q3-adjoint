#!/usr/bin/env python3
"""Koszul complex for F = tensor_j Sym^{k_j}(C^2_x + C^2_y) and m fundamental minus slots.
C^t = sum_{|T|=t} F (x) (C2_x)^{T^c} (x) (C2_y)^T, d moves one slot x->y using the abelian
gl(4) block Hom(C2_x,C2_y): Phi(e+ (x) P) = f- (x) f+ d_{e-}P - f+ (x) f- d_{e-}P, etc.
Reports SU(2)_y-invariant cohomology as SU(2)_x characters per degree (mod-p ranks).
usage: probe_su2_fm_nkoszul.py m k1 k2 ...   (e.g. probe_su2_fm_nkoszul.py 2 2)"""
import sys, itertools
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
m=int(sys.argv[1]); args=sys.argv[2:]
coef=None
if "--coef" in args:
    i=args.index("--coef"); coef=[int(c) for c in args[i+1].split(",")]; args=args[:i]
ks=[int(a) for a in args]
PR=2**31-1
J=len(ks)
if coef is None: coef=[1]*J
# variables: per factor j: e+ e- f+ f-  (index 4j+0..3); weights: e+:(1,0) e-:(-1,0) f+:(0,1) f-:(0,-1)
VW=[(1,0),(-1,0),(0,1),(0,-1)]
def monos_deg(k):
    for a in range(k+1):
        for b in range(k+1-a):
            for c in range(k+1-a-b):
                yield (a,b,c,k-a-b-c)
Fbasis=[tuple(itertools.chain(*combo)) for combo in itertools.product(*[list(monos_deg(k)) for k in ks])]
def wt_F(mono):
    x=y=0
    for j in range(J):
        a,b,c,d=mono[4*j:4*j+4]; x+=a-b; y+=c-d
    return (x,y)
# slot states: each slot is ('x',+1/-1) or ('y',+1/-1)
def slot_wt(s): return (s[1],0) if s[0]=='x' else (0,s[1])
def basis_deg(t):
    out=[]
    for T in itertools.combinations(range(m),t):
        choices=[[('y',1),('y',-1)] if r in T else [('x',1),('x',-1)] for r in range(m)]
        for sl in itertools.product(*choices):
            for mono in Fbasis:
                out.append((sl,mono))
    return out
def wt(el):
    sl,mono=el; x,y=wt_F(mono)
    for s in sl: dx,dy=slot_wt(s); x+=dx; y+=dy
    return (x,y)
def apply_d(el):
    """returns dict el'->coef"""
    sl,mono=el; res={}
    ypos=0
    for r in range(m):
        if sl[r][0]=='y': ypos+=1; continue
        sign=(-1)**ypos  # Koszul sign: number of odd slots before r
        eps=sl[r][1]
        # v=e+ contracts e- variables (coef +1); v=e- contracts e+ (coef -1)
        dvar=1 if eps==1 else 0; c0=1 if eps==1 else -1
        for j in range(J):
            p=mono[4*j+dvar]
            if p==0: continue
            base=list(mono); base[4*j+dvar]-=1
            # invariant f+ (x) f- - f- (x) f+ : slot gets first, polynomial multiplied by second
            for (slotsign,mulvar,cc) in ((1,3,1),(-1,2,-1)):
                nm=list(base); nm[4*j+mulvar]+=1
                nsl=list(sl); nsl[r]=('y',slotsign)
                key=(tuple(nsl),tuple(nm))
                res[key]=(res.get(key,0)+sign*c0*cc*p*coef[j])%PR
    return {k:v for k,v in res.items() if v}
def rank_mod(rows,ncols_index):
    # rows: list of dict col->val
    piv={}; r=0
    for row in rows:
        row=dict(row)
        while row:
            c=min(row); v=row[c]
            if c in piv:
                prow=piv[c]; f=v
                for cc,vv in prow.items():
                    row[cc]=(row.get(cc,0)-f*vv)%PR
                    if row[cc]==0: del row[cc]
            else:
                inv=pow(v,PR-2,PR)
                piv[c]={cc:(vv*inv)%PR for cc,vv in row.items()}
                r+=1; break
    return r
B=[basis_deg(t) for t in range(m+1)]
byw=[{} for t in range(m+1)]
for t in range(m+1):
    for el in B[t]: byw[t].setdefault(wt(el),[]).append(el)
# sanity d^2=0 on a sample
for t in range(m-1):
    for el in B[t][:200]:
        tot={}
        for k1,v1 in apply_d(el).items():
            for k2,v2 in apply_d(k1).items(): tot[k2]=(tot.get(k2,0)+v1*v2)%PR
        assert all(v==0 for v in tot.values()), "d^2 != 0"
rk={}
def rank_d(t,w):
    if (t,w) in rk: return rk[(t,w)]
    if t<0 or t>=m: rk[(t,w)]=0; return 0
    src=byw[t].get(w,[]); tgt=byw[t+1].get(w,[])
    idx={el:i for i,el in enumerate(tgt)}
    rows=[{idx[k]:v for k,v in apply_d(el).items()} for el in src]
    rk[(t,w)]=rank_mod(rows,len(tgt)); return rk[(t,w)]
def hdim(t,w):
    return len(byw[t].get(w,[]))-rank_d(t,w)-rank_d(t-1,w)
maxx=sum(ks)+m
print("F=Sym^%s m=%d"%(ks,m))
for t in range(m+1):
    chars={}
    for a in range(0,maxx+1):
        mult=hdim(t,(a,0))-hdim(t,(a+2,0))-hdim(t,(a,2))+hdim(t,(a+2,2))
        if mult: chars[a]=mult
    print("  H^%d  G_y-invariant part as G_x-char: %s"%(t,chars),flush=True)
