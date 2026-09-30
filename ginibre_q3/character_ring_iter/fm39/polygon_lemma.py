# Geometric lemma test: points g_0..g_M in the plane with (i) g_k ^ g_(k+1) >= 0 (ccw about O) and
# (ii) left turns (g_k - g_(k-1)) ^ (g_(k+1) - g_k) >= 0.  Is the signed area of every sub-polygon g_i..g_j >= 0?
import random, math
def wedge(p,q): return p[0]*q[1]-p[1]*q[0]
def area2(pts): return sum(wedge(pts[k],pts[(k+1)%len(pts)]) for k in range(len(pts)))
random.seed(3)
tot=bad=0; ex=None
for trial in range(20000):
    # random ccw curve: polar steps with angle increments in [0,pi), radii random -> keep only locally convex ones
    M=random.randint(3,14); th=random.uniform(0,6.3); pts=[]
    rad=random.uniform(.2,5)
    for k in range(M+1):
        pts.append((rad*math.cos(th),rad*math.sin(th)))
        th+=random.uniform(0,2.0); rad*=math.exp(random.uniform(-.8,.8))
    ok=all(wedge(pts[k],pts[k+1])>=0 for k in range(M)) and all(wedge((pts[k][0]-pts[k-1][0],pts[k][1]-pts[k-1][1]),(pts[k+1][0]-pts[k][0],pts[k+1][1]-pts[k][1]))>=0 for k in range(1,M))
    if not ok: continue
    tot+=1
    for i in range(M+1):
        for j in range(i+2,M+1):
            if area2(pts[i:j+1])< -1e-9:
                bad+=1; ex=ex or (pts[i:j+1],area2(pts[i:j+1])); break
        else: continue
        break
print('locally convex ccw curves tested',tot,'with a negative arc-chord area',bad)
if ex: print('example: total angle',round(sum(math.atan2(wedge(ex[0][k],ex[0][k+1]),ex[0][k][0]*ex[0][k+1][0]+ex[0][k][1]*ex[0][k+1][1]) for k in range(len(ex[0])-1)),3),'area2',ex[1])
print('--- split by angular sweep of the sub-arc ---')
def sweep(pts): return sum(math.atan2(wedge(pts[k],pts[k+1]),pts[k][0]*pts[k+1][0]+pts[k][1]*pts[k+1][1]) for k in range(len(pts)-1))
random.seed(4); stat={'sweep<=2pi':[0,0],'sweep>2pi':[0,0]}
for trial in range(60000):
    M=random.randint(3,14); th=random.uniform(0,6.3); pts=[]; rad=random.uniform(.2,5)
    for k in range(M+1):
        pts.append((rad*math.cos(th),rad*math.sin(th))); th+=random.uniform(0,2.0); rad*=math.exp(random.uniform(-.8,.8))
    ok=all(wedge(pts[k],pts[k+1])>=0 for k in range(M)) and all(wedge((pts[k][0]-pts[k-1][0],pts[k][1]-pts[k-1][1]),(pts[k+1][0]-pts[k][0],pts[k+1][1]-pts[k][1]))>=0 for k in range(1,M))
    if not ok: continue
    for i in range(M+1):
        for j in range(i+2,M+1):
            sub=pts[i:j+1]; key='sweep<=2pi' if sweep(sub)<=2*math.pi else 'sweep>2pi'
            stat[key][0]+=1; stat[key][1]+= area2(sub)< -1e-9
print(stat)
