import os, subprocess, sys, re
from collections import Counter

src=open('ginibre_q3/character_ring_iter/fm39/sec162_gp_single.cpp').read()
fd=os.memfd_create('sec163_witnesses',0)
os.set_inheritable(fd,True)
env=dict(os.environ)
env['TMPDIR']='/dev/shm'
b=subprocess.run(
    ['g++','-pipe','-std=c++17','-O3','-x','c++','-','-lgmpxx','-lgmp',
     '-o',f'/proc/self/fd/{fd}'],
    input=src.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,
    pass_fds=(fd,),env=env)
if b.returncode:
    print(b.stderr.decode(),file=sys.stderr)
    raise SystemExit(b.returncode)
os.fchmod(fd,0o700)

def add(B,n,m,s):
    B.extend([s*n]*m)

def profiles():
    B=[]
    for n in range(1,16,2):
        add(B,n,2,-1 if n==1 else 1)
    yield 'W128',16,B

    B=[-1]+[-2]*2
    for n in range(3,16,2):
        add(B,n,2,-1)
    yield 'W131',15,B

    yield 'W34',6,[-1]*4+[-2]+[-3]*8+[-4]

    B=[-1]*10+[-8]*4
    for n in range(3,32,2):
        add(B,n,31//n,-1)
    yield 'W427',31,B

    for sign,name in [(-1,'W540-'),(1,'W540-reflected')]:
        B=[]
        for n in range(1,36,2):
            if n==9:
                add(B,n,4,sign)
            else:
                add(B,n,35//n,-sign)
        yield name,36,B

    vals=[(-1,62),(-3,20),(-5,12),(-7,8),(9,7),(-11,5),(-13,4),
          (-15,4),(-17,3),(-19,2),(-21,2),(-23,2),(-25,2),(-27,2),
          (-29,2),(-31,2)]
    B=[x for x,m in vals for _ in range(m)]
    yield 'W869',57,B

    B=[-11]*5+[-54,-2]
    for n in range(1,54,2):
        if n!=11:
            add(B,n,54//n,1)
    yield 'W1283',55,B

def make_candidates(B):
    c=Counter(B)
    labs=sorted(c,key=lambda x:(abs(x),x))

    def valid(R):
        if not R:
            return False
        cc=Counter(R)
        return (sum(abs(x) for x in R)%2==0
                and all(cc[x]<=c[x] for x in cc))

    def rule(v):
        if c[v]>=2:
            return (v,v)
        same=[x for x in labs if x!=v and (abs(x)+abs(v))%2==0]
        ev=[x for x in labs if abs(x)%2==0]
        if same:
            return (v,max(same,key=lambda x:(abs(x),x)))
        if abs(v)%2==0:
            return (v,)
        if ev:
            return (min(ev,key=lambda x:(abs(x),x)),)
        return ()

    reps=[x for x in labs if c[x]>=2]
    pairs=[]
    for x in labs:
        for y in labs:
            if (abs(x)+abs(y))%2==0 and (x!=y or c[x]>=2):
                pairs.append((x,y))
    pairs.sort(key=lambda R:(abs(R[0]),R[0],abs(R[1]),R[1]))
    rw=max(labs,key=lambda x:(c[x]*abs(x),abs(x)))
    mx=max(labs,key=lambda x:(abs(x),x))
    most=max(labs,key=lambda x:(c[x],abs(x),x))
    same=[x for x in labs if x!=mx and (abs(x)+abs(mx))%2==0]
    ev=[x for x in labs if abs(x)%2==0]
    top=[]
    for x in sorted(B,key=lambda z:(abs(z),z),reverse=True):
        if not top:
            top=[x]
        elif (abs(x)+abs(top[0]))%2==0:
            top.append(x)
            break

    raw=[
        rule(rw), rule(mx),
        (min(reps,key=lambda x:(abs(x),x)),)*2 if reps else (),
        (max(reps,key=lambda x:(abs(x),x)),)*2 if reps else (),
        (rw,rw) if c[rw]>=2 else (),
        (most,most) if c[most]>=2 else (),
        pairs[0] if pairs else (), pairs[-1] if pairs else (),
        (mx,max(same,key=lambda x:(abs(x),x)))
            if abs(c[mx])==1 and same else (),
        (min(ev,key=lambda x:(abs(x),x)),) if ev else (),
        (max(ev,key=lambda x:(abs(x),x)),) if ev else (),
        tuple(top) if len(top)==2 else ()
    ]
    return [tuple(R) if valid(R) else None for R in raw]

names=['RuleW','RuleM','dup_small','dup_large','dup_heaviest',
       'dup_most_frequent','pair_smallest','pair_largest',
       'max_singleton_largest_same_parity','even_small','even_large',
       'top_two_compatible']

for name,p,B in profiles():
    cs=make_candidates(B)
    unique=[]
    for R in cs:
        if R is not None and R not in unique:
            unique.append(R)
    args=[str(p)]+list(map(str,B))
    for R in unique:
        args+=['--']+list(map(str,R))
    rr=subprocess.run([f'/proc/self/fd/{fd}']+args,pass_fds=(fd,),
                      env=env,capture_output=True,text=True,check=True)
    lines=rr.stdout.splitlines()
    g=int(lines[0].split('=')[1].strip())
    children={}
    for line,R in zip(lines[1:],unique):
        child=int(re.search(r'g_p\(B-R\) = (-?\d+)',line).group(1))
        children[R]=child
    mask=0
    print('WITNESS',name,'W',sum(map(abs,B)),'p',p,'factors',len(B),'g',g)
    for i,R in enumerate(cs):
        if R is not None:
            delta=g-children[R]
            if delta>=0:
                mask|=1<<i
            print(i,names[i],R,'delta',delta,'monotone',delta>=0)
    print('mask',mask)
