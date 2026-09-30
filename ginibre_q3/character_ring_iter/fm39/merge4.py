# Four-factor double merge: phi(h_A h_B h_C h_D) = sum_(d in CG(A,B), d' in CG(C,D)) phi(hatS_d hatS_d') - corr_pairing.
# If some pairing has corr <= 0, then (E) (T + W >= 0 for hatS hatS) gives phi >= 0.
exec(open('general_row_words.py').read().split("random.seed(")[0])
from collections import Counter
import itertools, random
def cg_(p,q): return range(abs(p-q),p+q+1,2)
random.seed(4); st=Counter(); ex=[]
for trial in range(1500):
    r=random.randint(2,5); e=2*r-4; a=random.randint(0,12)
    Ls=sorted([random.randint(3,9) for _ in range(4)],reverse=True)
    c=mkrow(a,e,[]); N=len(c)-1
    if (N+sum(Ls))%2: continue
    phi=phi_row(c,Ls,[],e%2)
    corrs=[]
    for (i,j),(k,l) in [((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]:
        main=sum(phi_row(c,[],[d,dp],e%2) for d in cg_(Ls[i],Ls[j]) for dp in cg_(Ls[k],Ls[l]))
        corrs.append(main-phi)
    st['words']+=1; st['phi<0']+= phi<0
    if min(corrs)<=0: st['some pairing corr<=0']+=1
    if min(corrs)==0: st['some pairing corr=0']+=1
    if min(corrs)>0 and len(ex)<4: ex.append((r,a,Ls,[float(x) for x in corrs],float(phi)))
print(dict(st)); print('all-positive-correction examples:',ex)
