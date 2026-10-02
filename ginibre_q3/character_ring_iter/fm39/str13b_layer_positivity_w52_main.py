# FM-STR13b (main agent): layer positivity for some pair on the full W <= 52 no-flip census (75,532 rows).
# Expected final line: TopPair ok 74788, fallback 744, no pair 0 (about 6 minutes).
import re, gzip, sys, time
import pathlib
_here=pathlib.Path(__file__).resolve().parent
src=open(_here/"str13_layer_positivity_main.py").read(); src=src[:src.index("LOG=pathlib")]
ns={}; exec(src,ns); layers,allpos,toppair=ns['layers'],ns['allpos'],ns['toppair']
LOG=_here/"sec166_census_w52_noflip.log.gz"
rows=[]
for line in gzip.decompress(open(LOG,'rb').read()).decode().splitlines():
    if not line.startswith("NOFLIP"): continue
    p=int(re.search(r"p=(-?\d+)",line).group(1)); B=[int(x) for x in re.search(r"B=([-\d ]+?)\s+phi",line).group(1).split()]
    rows.append(B+[p])
t0=time.time(); hb=time.time(); top_ok=0; fallback=0; none=[]
for k,L in enumerate(rows,1):
    i,j=toppair(L[:-1])
    if allpos(layers(L,i,j)): top_ok+=1
    elif any(allpos(layers(L,a,b)) for a in range(len(L)) for b in range(a+1,len(L))): fallback+=1
    else: none.append(L); print("NO PAIR",L,flush=True)
    if time.time()-hb>60: print(time.strftime('%H:%M:%S'),f"rows {k}/{len(rows)} TopPair ok {top_ok} fallback {fallback} none {len(none)}",flush=True); hb=time.time()
print(time.strftime('%H:%M:%S'),"FINAL rows",len(rows),"TopPair ok",top_ok,"fallback",fallback,"no pair",len(none),f"{time.time()-t0:.0f}s",flush=True)
