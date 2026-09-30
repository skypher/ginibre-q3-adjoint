# Which row class does (E) need?  Reciprocal rows (e even -> c_(N-k) = c_k): real-rooted |t|>=2, unit-circle |t|<2,
# complex quadruples with log-concave weight (|al|-2 >= be) and without.
exec(open('recip_test.py').read().split("def tools")[0].split("import random")[1].replace("from fractions import Fraction as Fr\nfrom math import comb\n",""))
import random
from fractions import Fraction as Fr
from collections import Counter
def Echeck(c):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); B=lambda k: C_(k-1)+C_(k+1)
    bad=0; tot=0
    for j in range((N+1)//2,N+1):
        for i in range(j+1,N+2):
            W=B(i)*C_(j)-C_(i)*B(j); tot+=1; bad+= D(j)-D(i)<abs(W)
    return tot,bad
random.seed(3); st=Counter()
for trial in range(400):
    e=random.choice([0,2,4]); a=random.randint(0,6)
    for fam in ['real |t|>=2','unit circle |t|<2','quartic log-concave','quartic not log-concave']:
        c=mkrow(a,e,[])
        if fam=='real |t|>=2':
            for _ in range(random.randint(1,2)): c=mulrow(c,[1,random.choice([-1,1])*(2+Fr(random.randint(0,20),5)),1])
        elif fam=='unit circle |t|<2':
            for _ in range(random.randint(1,2)): c=mulrow(c,[1,Fr(random.randint(-19,19),10),1])
        else:
            be=Fr(random.randint(1,30),10)
            al=(2+be+Fr(random.randint(0,20),10))*random.choice([-1,1]) if fam=='quartic log-concave' else Fr(random.randint(-15,15),10)
            if fam=='quartic not log-concave' and abs(al)-2>=be: al=Fr(0)
            c=mulrow(c,[1,2*al,2+al*al+be*be,2*al,1])
        tot,bad=Echeck(c); st[(fam,'pairs')]+=tot; st[(fam,'violations')]+=bad; st[(fam,'rows with violation')]+= bad>0
for k in sorted(st): print(k,st[k])
