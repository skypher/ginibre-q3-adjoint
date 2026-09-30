# Main agent: on 0 <= a,e <= 70, every consumer row of the q=2 plus slack has a PD recurrence quadratic or W >= 0.  Run from the repo root.
import sys
sys.argv=['x','--amax','1']
exec(open('ginibre_q3/character_ring_iter/fm39/sec127_q2plus_audit_repro.py').read().split("if __name__")[0])
amax=int(sys.argv_real[1]) if False else 70
nonpd=0; nonpd_Wneg=0; tot=0; ex=[]
for a in range(amax+1):
    for e in range(amax+1):
        n=a+e
        if n<2: continue
        c=row(a,e)
        for j in range((n+1)//2, n):
            A,B,C=quadratic(n,a-e,j); tot+=1
            pd = A>0 and B*B-4*A*C<0
            if not pd:
                nonpd+=1
                val,W=slack_and_W(c,n,j)
                if W<0: nonpd_Wneg+=1; ex.append((a,e,j,W,val))
                elif len(ex)<3: pass
print("rows",tot,"non-PD",nonpd,"non-PD with W<0",nonpd_Wneg,ex[:5])
