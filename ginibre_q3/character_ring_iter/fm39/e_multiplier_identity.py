# Guess: (E)_(q,-/+) at (j, i=j+q+1) = sum_(i0=0..floor((q-1)/2)) delta_(j+q)[(z^i0 -/+ z^(q-i0)) c]  (+ middle term for even q?)
# Search the exact index/multiplier pattern on random sequences.
import random, itertools
def rowmul(c,R):
    o=[0]*(len(c)+len(R)-1)
    for i,x in enumerate(c):
        for j,y in enumerate(R): o[i+j]+=x*y
    return o
def D(c,k):
    C_=lambda t: c[t] if 0<=t<len(c) else 0
    return C_(k)**2-C_(k-1)*C_(k+1)
def dl(c,k): return D(c,k)-D(c,k+1)
random.seed(1)
for q in (1,2,3,4,5):
    for sg in (1,-1):
        found=None
        for shift in range(0,q+3):
            ok=True
            for _ in range(30):
                c=[0]*3+[random.randint(-9,9) for _ in range(12)]+[0]*3
                j=8; i=j+q+1; C_=lambda t: c[t] if 0<=t<len(c) else 0; B=lambda k: C_(k-1)+C_(k+1)
                W=B(i)*C_(j)-C_(i)*B(j); tgt=D(c,j)-D(c,i)-sg*W
                tot=0
                for i0 in range(0,q//2+1):
                    if 2*i0>q: break
                    R=[0]*(q+1); R[i0]+=1; R[q-i0]+= -sg
                    if 2*i0==q:   # middle term: (1 - sg) z^(q/2): zero for sg=+1, 2 z^(q/2) for sg=-1
                        R=[0]*(q+1); R[i0]= (1-sg)
                        w=0.5  # half weight guess
                    else: w=1
                    tot+=w*dl(rowmul(c,R),j+shift)
                if tot!=tgt: ok=False; break
            if ok: found=shift; break
        print('q=%d sign=%+d: identity with shift'%(q,sg), found)
