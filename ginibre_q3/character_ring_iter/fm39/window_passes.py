# Theorem WL test: for windows with sweep > 2 pi and T > 0, check whether the first pass lies outside g_end at direction(g_end)
# or the last pass lies outside g_x at direction(g_x).  If one holds, the triangle O,g_x,g_end sits in that pass's fan.
import math
exec(open('split3.py').read())
def ang(p,q): return math.atan2(p[0]*q[1]-p[1]*q[0], p[0]*q[0]+p[1]*q[1])
def ray_radius(p,q,phi):   # radius where segment p->q meets the ray at absolute direction phi
    d=(math.cos(phi),math.sin(phi)); ex=(q[0]-p[0],q[1]-p[1])
    den=d[0]*ex[1]-d[1]*ex[0]
    if abs(den)<1e-300: return math.hypot(*p)
    t=(p[0]*ex[1]-p[1]*ex[0])/den
    return t
def pass_radius(g,beta,phi):  # walk from g[0]; first edge where cumulative angle reaches beta; radius at direction phi
    cum=0.0
    for k in range(len(g)-1):
        dth=ang(g[k],g[k+1])
        if cum+dth>=beta-1e-15: return ray_radius(g[k],g[k+1],phi)
        cum+=dth
    return None
stat={'long&T>0':0,'first ok':0,'last ok':0,'neither':0}; ex=None
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=list(cv(a,e)); C_=lambda k: c[k] if 0<=k<=N else 0
        for x in range(0,N+1):
            for Cw in range(1,N-x+1):
                g=[(float(C_(k)),float(C_(k-1))) for k in range(x,x+Cw+2)]
                sw=sum(ang(g[k],g[k+1]) for k in range(len(g)-1))
                T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
                if sw<=2*math.pi or T<=0: continue
                stat['long&T>0']+=1
                beta=sw%(2*math.pi); thx=math.atan2(g[0][1],g[0][0]); the=math.atan2(g[-1][1],g[-1][0])
                r1=pass_radius(g,beta,the); rL=pass_radius(g[::-1],beta,thx)   # reversed walk is clockwise; angles negative
                if rL is None or True:
                    # reversed curve turns clockwise: redo with negated angles
                    gr=[(p[0],-p[1]) for p in g[::-1]]; rL=pass_radius(gr,beta,-thx)
                f1 = r1 is not None and r1>=math.hypot(*g[-1])*(1-1e-12)
                fL = rL is not None and rL>=math.hypot(*g[0])*(1-1e-12)
                stat['first ok']+=f1; stat['last ok']+=fL
                if not(f1 or fL): stat['neither']+=1; ex=ex or (r,a,x,Cw,sw)
print(stat,'first "neither" case:',ex)
