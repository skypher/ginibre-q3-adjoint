import argparse,datetime,itertools,math
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache

ap=argparse.ArgumentParser(
    description='Exact FM-STR4f verifier for interior-prefix identities and separators.')
ap.parse_args()

def log(*x):
    print(datetime.datetime.now(datetime.timezone.utc).strftime('%H:%M:%S UTC'),*x,flush=True)

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def char(w):
    d={(0,0):1}
    for n,e in w:
        z=defaultdict(int)
        for (r,s),v in d.items():
            for q in cg(r,n):
                z[q,s]+=v
            for q in cg(s,n):
                z[r,q]+=e*v
        d={k:v for k,v in z.items() if v}
    return d

def phi(w):
    return char(tuple(w)).get((0,0),0)

def top_pair(w):
    return max(
        (w[i][0]+w[j][0],max(w[i][0],w[j][0]),-i,-j,i,j)
        for i in range(len(w)) for j in range(i+1,len(w))
        if (w[i][0]+w[j][0])%2==0
    )[-2:]

def split_rest(w,i,j):
    C=tuple(x for k,x in enumerate(w) if k not in (i,j))
    order=sorted(range(len(C)),key=lambda k:-C[k][0])
    A=[]; B=[]; wa=wb=0
    for k in order:
        if wa<=wb:
            A.append(k); wa+=C[k][0]
        else:
            B.append(k); wb+=C[k][0]
    if wa>wb:
        A,B=B,A
    return C,tuple(A),tuple(B)

@lru_cache(None)
def prefix_layers(w,i,j):
    C,Ai,Bi=split_rest(w,i,j)
    A=tuple(C[k] for k in Ai)
    B=tuple(C[k] for k in Bi)+(w[i],w[j])
    fa,fb=char(A),char(B)
    lay=[0]*(sum(n for n,e in A)+1)
    for (r,s),v in fa.items():
        lay[r+s]+=v*fb.get((r,s),0)
    for t in range(1,len(lay)):
        lay[t]+=lay[t-1]
    return tuple(lay)

def pi(w,i,j,T):
    z=prefix_layers(tuple(w),i,j)
    return z[min(T,len(z)-1)]

def even_signs(n):
    return [s for s in itertools.product((1,-1),repeat=n) if math.prod(s)==1]

def word(L,s):
    return tuple(zip(L,s))

def fuse(L,s,i,j,c):
    return tuple((L[k],s[k]) for k in range(len(L)) if k not in (i,j))+((c,s[i]*s[j]),)

def child_top_pi(w,T):
    return pi(w,*top_pair(w),T)

def mval(L):
    d={0:1}
    for n in L:
        z=defaultdict(int)
        for r,v in d.items():
            for q in cg(r,n):
                z[q]+=v
        d={k:v for k,v in z.items() if v}
    return d.get(0,0)

def noflip(L):
    out=[]
    for s in even_signs(len(L)):
        v=phi(word(L,s)); ds=[]
        for i,j in itertools.combinations(range(len(L)),2):
            t=list(s); t[i]*=-1; t[j]*=-1
            ds.append(v-phi(word(L,t)))
        if all(d<0 for d in ds):
            out.append((s,v,max(ds)))
    return out

def candidate_columns(L):
    S=even_signs(len(L)); C=[]; negative=0
    # Fusion children: TopPair cut, all channels, all finite heights.
    for i,j in itertools.combinations(range(len(L)),2):
        for c in cg(L[i],L[j]):
            ws=[fuse(L,s,i,j,c) for s in S]
            vv=[prefix_layers(w,*top_pair(w)) for w in ws]
            for T in range(max(map(len,vv))):
                col=[v[min(T,len(v)-1)] for v in vv]
                if any(col):
                    C.append(col); negative+=min(col)<0
    # Delete equal-sign pairs; test every child cut and height.
    for i,j in itertools.combinations(range(len(L)),2):
        ws=[
            (tuple((L[k],s[k]) for k in range(len(L)) if k not in (i,j))
             if s[i]*s[j]==1 else None)
            for s in S
        ]
        for u,v in itertools.combinations(range(len(L)-2),2):
            ids=[r for r,w in enumerate(ws) if w is not None]
            vv=[prefix_layers(ws[r],u,v) for r in ids]
            if not vv:
                continue
            for T in range(max(map(len,vv))):
                col=[0]*len(S)
                for q,r in enumerate(ids):
                    col[r]=vv[q][min(T,len(vv[q])-1)]
                if any(col):
                    C.append(col)
                    negative+=min(x for x in col if x!=0)<0
    # If m(L)>0, m(empty)m(L) spans the nonnegative constant ray.
    if mval(L)>0:
        C.append([1]*len(S))
    return S,C,negative

Y_PRIOR={
0:Q(-8164,182391),2:Q(-3195,549569),8:Q(-20673,997262),
10:Q(11083,336197),11:Q(-2725,111192),12:Q(23493,947032),
14:Q(1975,234776),18:Q(12803,303585),21:Q(26933,543901),
23:Q(-15007,609018),25:Q(3197,196707),26:Q(27109,464807),
27:Q(-24023,387090),29:Q(-81969,999596),30:Q(-30081,879023),
33:Q(6260,306561),34:Q(14597,325627),36:Q(18652,606331),
37:Q(-20255,894823),38:Q(-40907,972609),40:Q(8983,813068),
43:Q(-3053,800052),45:Q(-10668,618607),46:Q(9271,839151),
49:Q(-34607,885173),51:Q(-1486,61861),53:Q(-27827,682774),
55:Q(1894,78971),57:Q(19670,800447),60:Q(-1354,283449),
63:Q(50779,471301)
}

Y_PLAIN8={
0:Q(-121454,882437),1:Q(1528,488793),5:Q(25143,602846),
6:Q(41922,926053),7:Q(7919,741872),8:Q(12353,654692),
10:Q(85493,924934),11:Q(9,3590),12:Q(49942,898527),
13:Q(-9177,259688),15:Q(-48795,960851),16:Q(-13131,584852),
17:Q(27533,518082),19:Q(-1449,297604),21:Q(-5983,198031),
22:Q(-8152,987611),27:Q(-8781,498706),28:Q(5157,310699),
30:Q(-15151,836907),31:Q(2234,960973),33:Q(5015,578808),
37:Q(-9304,931715),39:Q(-8075,358414),46:Q(-39471,997984),
47:Q(25253,999657),48:Q(27567,533857),51:Q(-2647,668412),
52:Q(-5239,158488),57:Q(-2104,76945),60:Q(-1486,523219),
62:Q(8797,519183),63:Q(49717,547920)
}

def verify_sep(name,L,Y):
    S,C,negative=candidate_columns(L)
    P=[phi(word(L,s)) for s in S]
    score=sum(Y.get(r,0)*P[r] for r in range(len(S)))
    dots=[sum(Y.get(r,0)*col[r] for r in range(len(S))) for col in C]
    assert score<0 and all(x>=0 for x in dots) and negative==0
    log('EXACT_SEPARATOR',name,'rows',len(S),'columns',len(C),
        'negative-child-columns',negative,'y-support',len(Y),
        'target-dot',score,'min-column-dot',min(dots))

log('START exact re-expansion')
prior=(1,4,5,6,6,8,9,15)
census=(1,3,3,4,5,5,6,7)
verify_sep('prior obstruction; fused child TopPair, deletions all cuts',
           prior,Y_PRIOR)
verify_sep('plain run (1,...,8); fused child TopPair, deletions all cuts',
           tuple(range(1,9)),Y_PLAIN8)

for name,L,v in [('prior',prior,3532),('census',census,1100)]:
    nf=noflip(L); M=mval(L); gamma=Q(v,M)
    assert len(nf)==2 and M>0 and all(x==v for s,x,d in nf)
    assert all(gamma*M==x for s,x,d in nf)
    log('NOFLIP_CONSTANT',name,
        [(''.join('+' if e>0 else '-' for e in s),x,d) for s,x,d in nf],
        'm(empty)m(L)',M,'gamma',gamma)

def run_identity(k,p,support,count,phis):
    L=tuple(range(1,k+1))+(p,)
    W=k*(k+1)//2; delta=(W-p)//2
    assert p>k and p>=6 and (W-p)%2==0 and delta>=8
    expected=-1 if k%2 else 1
    rows=noflip(L)
    assert len(rows)==count and [v for s,v,d in rows]==phis
    assert any(tuple(s)==tuple([-1]*k+[expected]) for s,v,d in rows)
    for s,v,d in rows:
        rhs=Q(0)
        for term in support:
            if term[0]=='F':
                _,a,i,j,c,T=term
                rhs+=a*child_top_pi(fuse(L,s,i,j,c),T)
            else:
                _,a,i,j,u,vv,T=term
                if s[i]*s[j]==1:
                    child=tuple((L[q],s[q]) for q in range(len(L))
                                if q not in (i,j))
                    rhs+=a*pi(child,u,vv,T)
        assert rhs==v,(k,p,s,rhs,v)
    log('RESIDUAL_RUN',k,p,'delta',delta,'noflip_rows',count,
        'Phi',phis,'support',
        [(z[0],str(z[1]),z[2:]) for z in support],'PASS')

RUNS=[
(7,8,[('F',Q(239,12),0,3,5,0),('D',Q(490,11),0,2,0,1,3)],
 4,[956,980,956,980]),
(7,10,[('F',Q(457,238),0,1,3,8),('D',Q(1466,357),0,2,0,1,1)],
 4,[914,954,914,954]),
(7,12,[('F',Q(1311,476),2,6,10,9),('F',Q(199,119),5,6,1,6)],
 4,[806,830,806,830]),
(8,10,[('F',Q(3217,38),0,1,1,0),('F',Q(3105,38),0,5,7,0)],
 4,[6210,6434,6210,6434]),
(8,12,[('F',Q(2953,36),0,1,3,0),('F',Q(1015,12),0,3,5,0)],
 4,[5906,6090,5906,6090]),
(8,14,[('F',Q(97),0,1,1,0),('F',Q(1321,14),0,5,7,0)],
 4,[5284,5432,5284,5432]),
(8,16,[('F',Q(568,5),0,1,1,0),('F',Q(1113,10),0,1,3,0),
       ('F',Q(371),5,6,5,0)],
 6,[4452,4452,4544,4452,4452,4544]),
(8,18,[('F',Q(1789,12),0,1,1,0),('F',Q(1745,12),1,4,7,0)],
 4,[3490,3578,3490,3578]),
(8,20,[('F',Q(1315,6),0,1,1,0),('F',Q(1291,6),1,4,7,0)],
 4,[2582,2630,2582,2630])
]
for args in RUNS:
    run_identity(*args)
log('ALL_EXACT_CHECKS_PASS')
