from random import Random
rng=Random(169108)

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

def phi(word):
    d={(0,0):1}
    for e,n in word:
        z={}
        for (a,b),v in d.items():
            for x in cg(a,n):
                z[x,b]=z.get((x,b),0)+v
            for y in cg(b,n):
                z[a,y]=z.get((a,y),0)+e*v
        d=z
    return d.get((0,0),0)

for _ in range(250):
    C=[(rng.choice((-1,1)),rng.randrange(1,10))
       for _ in range(rng.randrange(8))]
    n=rng.randrange(1,9)
    assert phi(C+[(1,n),(-1,n)]) == sum(
        phi(C+[(-1,2*k)]) for k in range(1,n+1))

for _ in range(250):
    C=[(rng.choice((-1,1)),2*rng.randrange(1,8))
       for _ in range(rng.randrange(7))]
    a=2*rng.randrange(8)+1
    b=2*rng.randrange(8)+1
    e=rng.choice((-1,1))
    f=rng.choice((-1,1))
    assert phi(C+[(e,a),(f,b)]) == sum(
        phi(C+[(e*f,j)]) for j in cg(a,b))

for _ in range(300):
    w=[(rng.choice((-1,1)),rng.randrange(1,12))
       for _ in range(rng.randrange(1,9))]
    if sum(n for e,n in w)%2 or sum(e<0 for e,n in w)%2:
        assert phi(w)==0

print("pair=250 two_odd=250 parity=300 PASS")
