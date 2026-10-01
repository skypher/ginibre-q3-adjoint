import re, subprocess, sys, time
from fractions import Fraction
from concurrent.futures import ThreadPoolExecutor
GP='../dtest/gp_single'
rows=[]
for l in open(sys.argv[1]):
    if not l.startswith('NOFLIP'): continue
    m=re.match(r'NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+).*TopPair\((-?\d+),(-?\d+)\)',l)
    rows.append((int(m.group(1)),abs(int(m.group(2))),m.group(3).split(),m.group(5),m.group(6)))
def run(r):
    W,p,B,a,b=r
    out=subprocess.run([GP,str(p)]+B+['--',a,b],capture_output=True,text=True).stdout
    g=int(re.search(r'g_p\(B\) = (-?\d+)',out).group(1)); h=int(re.search(r'g_p\(B-R\) = (-?\d+)',out).group(1))
    return (Fraction(h,g),W,p,B)
t0=time.time(); best={}
with ThreadPoolExecutor(16) as ex:
    for k,(fr,W,p,B) in enumerate(ex.map(run,rows,chunksize=64)):
        if W not in best or fr>best[W][0]: best[W]=(fr,p,B)
        if k%5000==0: print(time.strftime('%H:%M:%S'),'done',k,'/',len(rows),flush=True)
for W in sorted(best): print('W',W,'max child/parent %.4f'%float(best[W][0]),'p',best[W][1],'B',' '.join(best[W][2]))
