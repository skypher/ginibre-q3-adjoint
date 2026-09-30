# Independent check of the parity transformation (1) and the three-factor reduction (4), with the main agent's phi.
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
import itertools
bad=0; n=0
for r in range(1,5):
    for m in range(1,min(2*r,4)+1):
        for kap in itertools.combinations_with_replacement(range(1,6),m):
            for a in range(0,4):
                lhs_w=ONE
                for k in kap: lhs_w=mul(lhs_w,h(k))
                L=phi(r,lhs_w,a)
                odd=[k for k in kap if k%2]; even=[k for k in kap if k%2==0]
                t=len(odd); e=2*r-m
                if (a+t)%2:
                    n+=1; bad+= (L!=0); continue
                R=(a+t)//2
                w=ONE
                for k in odd: w=mul(w,h(k))
                for k in even: w=mul(w,shat(k+1))
                if R>=1: Rv=phi(R,w,e)
                else: Rv=E(mul(w,dict(spow(e))))/2
                n+=1
                if L!=Rv: bad+=1; print('MISMATCH',r,kap,a,L,Rv) if bad<4 else None
print('parity transformation checks:',n,'mismatches',bad)
# three-factor sector: phi_r(h_u h_v h_w) >= 0 for r<=8, u,v,w<=7
neg=0; cnt=0
for r in range(1,9):
    for u,v,w in itertools.combinations_with_replacement(range(1,8),3):
        val=phi(r,mul(mul(h(u),h(v)),h(w)),0); cnt+=1
        if val<0: neg+=1
print('phi_r(h_u h_v h_w), r<=8, parts<=7:',cnt,'values, negatives',neg)
