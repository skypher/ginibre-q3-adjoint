from datetime import datetime, timezone
from fractions import Fraction
def log(*x):
    print(datetime.now(timezone.utc).isoformat(timespec="seconds"), *x, flush=True)
def mul(x,y):
    z={}
    for a,c in x.items():
        for b,d in y.items():
            k=tuple(i+j for i,j in zip(a,b))
            z[k]=z.get(k,Fraction())+c*d
    return {k:v for k,v in z.items() if v}
def sub(x,y):
    z=dict(x)
    for k,v in y.items(): z[k]=z.get(k,Fraction())-v
    return {k:v for k,v in z.items() if v}
def ent(c,s,n):
    return {tuple(int(i in s) for i in range(n)):Fraction(c)}
def det2(m):
    return sub(mul(m[0][0],m[1][1]),mul(m[0][1],m[1][0]))
log("START exact polynomial checks")
S={0,8,9}; Sc=set(range(12))-S
for a in range(4,65):
    h=Fraction((3*a+6)*(3*a+7),(a+1)*(2*a+6))
    assert h>0
    m=[[ent(h,S,12),ent(h,Sc,12)],
       [ent(h,set(),12),ent(h,set(),12)]]
    d=det2(m)
    assert len(d)==2
    assert set(d)=={tuple(int(i in S) for i in range(12)),
                    tuple(int(i in Sc) for i in range(12))}
    assert sorted(abs(c) for c in d.values())==[h*h,h*h]
log("PASS Lambda_a determinant factor, a=4..64")
m=[[ent(5,S,12),ent(5,Sc,12)],
   [ent(Fraction(5,2),S,12),ent(Fraction(5,2),Sc,12)]]
assert all(m[i][j] for i in range(2) for j in range(2))
assert det2(m)=={}
log("PASS K2,2 support with identically zero determinant")
for q in range(1,11):
    odd=[s for s in range(1<<q) if s.bit_count()%2]
    even=[t for t in range(1<<q) if t.bit_count()%2==0]
    assert len(odd)==len(even)==1<<(q-1)
    a=[[int((s^t)==1) for s in odd] for t in even]
    assert all(sum(row)==1 for row in a)
    assert all(sum(a[i][j] for i in range(len(even)))==1
               for j in range(len(odd)))
log("PASS cube leading permutation for q=1..10")
log("PASS all checks")
