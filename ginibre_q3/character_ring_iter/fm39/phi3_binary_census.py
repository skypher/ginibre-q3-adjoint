# Wider census: three-factor words (both branches), binary-form definiteness.  Heartbeat.
import time, sys
exec(open('phi3_binary_form.py').read().split("st=Counter(); chk=0")[0])
R2,A2=int(sys.argv[1]),int(sys.argv[2])
st=Counter(); t0=time.time(); last=t0; nondef=[]
for r in range(2,R2+1):
    e=2*r-3
    for a in range(0,A2+1):
        N=a+e
        if abs(a-e)<=1 or a==0: continue
        for w in range(2,min(N,14)+1):
            for v in range(w,N+3):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2
                    if not(al>=0): continue
                    br='G0' if al+B>=N+1 else 'gamma<=N'
                    Aq,Bq,Cq=phi_form(a,e,A,B,C); st[br]+=1
                    if Bq*Bq-4*Aq*Cq<0 and Aq>0: st[br+' definite']+=1
                    elif len(nondef)<8: nondef.append((r,a,u,v,w,br))
                    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat r=%d a=%d'%(r,a),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed',round(time.time()-t0,1)); print('non-definite examples:',nondef)
