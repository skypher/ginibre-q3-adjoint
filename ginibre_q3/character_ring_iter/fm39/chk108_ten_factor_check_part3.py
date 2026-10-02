from math import comb
N=10
M=1<<N

def mult0(a,mask):
    ids=[i for i in range(N) if mask>>i&1]
    r=len(ids)
    if not r:return 1
    w=sum(a[i] for i in ids)
    if r==1 or w%2:return 0
    z=0
    for v in range(1<<r):
        top=w//2+r-2-sum(a[ids[j]]+1 for j in range(r) if v>>j&1)
        z+=(-1 if v.bit_count()&1 else 1)*(comb(top,r-2) if top>=r-2 else 0)
    return z

def values(a):
    m=[mult0(a,s) for s in range(M)]
    f=[m[s]*m[M-1-s] for s in range(M)]
    h=1
    while h<M:
        for b in range(0,M,2*h):
            for j in range(b,b+h):
                x,y=f[j],f[j+h]
                f[j],f[j+h]=x+y,x-y
        h*=2
    return f

def valid(a):
    D=sum(a[:6])
    q,r,s,p=a[6:]
    v=D+q+r+s-p
    return (a[:6]==tuple(sorted(a[:6])) and a[5]<=8 and q>=a[5]
        and q<=r<=s<=p and p>=max(s,2*q-D,6)
        and p<=min(q+D,D+q+r-s) and v%2==0 and v//2>=8
        and s<=v//2
        and sum(x>=3 for x in a[:6])+(q>=3)+(r>=3)+(s>=3)>=2)

def signs(a):
    groups=[]
    for i,n in enumerate(a):
        if not groups or a[groups[-1][0]]!=n:groups.append([i])
        else:groups[-1].append(i)
    for u in range(1<<len(groups)):
        m=sum(1<<i for j,g in enumerate(groups) if u>>j&1 for i in g)
        if m.bit_count()%2==0:yield m

profiles={
    "max_p144_q96_M48_delta144":(8,8,8,8,8,8,96,144,144,144),
    "min_p6":(1,1,1,1,1,4,4,5,6,6),
    "q_2D_M_D":(1,1,1,1,1,8,26,39,39,39),
    "delta_s_q_c6":(1,1,1,1,1,8,8,8,8,21),
    "p_2q_minus_D":(1,1,1,1,1,8,18,23,23,23),
    "M_D":(1,1,1,1,1,8,18,31,31,31),
}
for name,a in profiles.items():
    assert valid(a),(name,a)
    t=values(a)
    ss=list(signs(a))
    low=min((t[m],m) for m in ss)
    D=sum(a[:6]);q,r,s,p=a[6:]
    print(name,"D",D,"M",p-q,"delta",(D+q+r+s-p)//2,
          "signings",len(ss),"min",low[0],"mask",low[1])
