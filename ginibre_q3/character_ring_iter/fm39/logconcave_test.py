# Complex t = al + i be in the factor (1+tz+z^2)(1+conj(t)z+z^2): weight K(x,2) gets ((x+al)^2+be^2);
# log-concave on [-2,2] iff |al|-2 >= be.  Do three-factor failures occur only when log-concavity fails?
exec(open('recip_bounds.py').read().split("for trial in range(600):")[0])
import random
random.seed(17)
res={'log-concave weight':[0,0,0],'weight not log-concave':[0,0,0]}
for trial in range(1500):
    e=random.choice([1,3,5]); a=random.randint(0,5)
    al=Fr(random.choice([-1,1])*random.randint(0,60),10); be=Fr(random.randint(1,40),10)
    quart=[1,2*al,2+al*al+be*be,2*al,1]
    c=mulrow(mkrow(a,e,[]),quart)
    key='log-concave weight' if abs(al)-2>=be else 'weight not log-concave'
    N=len(c)-1; rowneg=False
    for Cw in range(2,N+2):
        for B in range(Cw,N+3):
            for A in range(B,B+Cw+1):
                if (N+A-B-Cw)%2: continue
                v=Phi3(c,(N+A-B-Cw)//2,B,Cw); res[key][0]+=1
                if v<0: res[key][1]+=1; rowneg=True
    res[key][2]+=rowneg
for k,v in res.items(): print(k,': tests',v[0],'failures',v[1],'rows with failure',v[2])
