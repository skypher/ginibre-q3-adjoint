exec(open('recip_test.py').read().split("random.seed(11)")[0])
import random
random.seed(12)
stats={}
def rec(k,bad,ex=None):
    s=stats.setdefault(k,[0,0,0,None]); s[0]+=1
    if bad: s[1]+=1; s[3]=s[3] or ex
def three_all(c,tag,rowtag):
    N=len(c)-1; anyneg=False
    for Cw in range(2,N+2):
        for B in range(Cw,N+3):
            for A in range(B,B+Cw+1):
                if (N+A-B-Cw)%2: continue
                al=(N+A-B-Cw)//2; v=Phi3(c,al,B,Cw)
                rec(tag,v<0,(rowtag,A,B,Cw,str(v))); anyneg|= v<0
    s=stats[tag]; s[2]+=anyneg
for trial in range(600):
    # (i) non-reciprocal real-rooted: (1-z)^e times random real linear factors
    e=random.choice([1,3,5]); roots=[Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(1,5))]
    c=mkrow(0,e,[])
    for rho in roots: c=mulrow(c,[1,rho])
    three_all(c,'three-factor, NON-reciprocal real-rooted',(e,[str(x) for x in roots]))
    # (ii) reciprocal but unit-circle roots: |t|<2
    a=random.randint(0,5); ts=[Fr(random.randint(-19,19),10) for _ in range(random.randint(1,2))]
    c=mkrow(a,e,ts)
    three_all(c,'three-factor, reciprocal with unit-circle roots (|t|<2)',(a,e,[str(t) for t in ts]))
for k,v in stats.items(): print(k,': tests',v[0],'failures',v[1],'rows with a failure',v[2],'of 600; first',v[3])
