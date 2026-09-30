# LP: (E)_(q,sign) at (j, j+q+1) as a nonnegative combination of
#   (a) OL differences delta_(j+t)[R P] of base-family rows (2t > deg R), and
#   (b) proved (E) forms of lower q' (q' = 1, or q' = 2 minus sign) of base-family rows R P at (j+t, j+t+q'+1), valid when 2t >= deg R,
# as an identity of quadratic forms in a free sequence c.
import sys
import numpy as np
from scipy.optimize import linprog
exec(open('e_ol_lp.py').read().split("q=int(sys.argv[1])")[0])
def Eform(R,q,sg,j):
    i=j+q+1; c=lambda m: lin(R,m)
    Bi=addq(c(i-1),c(i+1)); Bj=addq(c(j-1),c(j+1))
    W=addq(prod(Bi,c(j),1),prod(c(i),Bj,-1))
    return addq(Dform(R,j),{k:-v for k,v in Dform(R,i).items()},{k:-sg*v for k,v in W.items()})
q=int(sys.argv[1]); dmax=int(sys.argv[2]); tmax=q+dmax+2
for sg in (1,-1):
    cols=[]; names=[]
    for d in range(0,dmax+1):
        for al in range(d+1):
            be=d-al; R=Rpoly(al,be)
            for t in range(-2, tmax+1):
                if 2*t>d: cols.append(deltaform(R,t)); names.append(('OL',al,be,t))
                if 2*t>=d:
                    for (qq,s2) in [(1,1),(1,-1),(2,1)]:
                        if qq<q: cols.append(Eform(R,qq,s2,t)); names.append(('E',qq,s2,al,be,t))
    tg=Etarget(q,sg)
    keys=sorted(set(k for cq in cols for k in cq)|set(tg))
    A=np.array([[float(cq.get(k,0)) for cq in cols] for k in keys]); b=np.array([float(tg.get(k,0)) for k in keys])
    res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
    print('q=%d sign=%+d dmax=%d cols=%d: %s'%(q,sg,dmax,len(cols),'FEASIBLE' if res.status==0 else 'infeasible'))
    if res.status==0: print('   ',[(n,round(w,4)) for w,n in zip(res.x,names) if w>1e-9][:12])
