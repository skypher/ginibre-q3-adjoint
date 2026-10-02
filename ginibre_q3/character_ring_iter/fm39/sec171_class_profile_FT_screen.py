import os, sys, re, time, random, subprocess
from fractions import Fraction

if any(x in ("-h", "--help") for x in sys.argv[1:]):
    print("Runs the FM-SEC171 97-profile exact class-pattern scan from the repo root.")
    raise SystemExit(0)

env = dict(os.environ, TMPDIR="/dev/shm", OMP_NUM_THREADS="4")

def compile_mem(name, path, openmp=True):
    fd = os.memfd_create(name, 0)
    os.set_inheritable(fd, True)
    flags = ["g++", "-pipe", "-std=c++17", "-O3"]
    if openmp:
        flags += ["-fopenmp"]
    flags += ["-x", "c++", "-o", f"/proc/self/fd/{fd}", "-",
              "-lgmpxx", "-lgmp"]
    with open(path, "rb") as f:
        src = f.read()
    q = subprocess.run(flags, input=src, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, pass_fds=(fd,), env=env)
    if q.returncode:
        raise RuntimeError(q.stderr.decode())
    os.fchmod(fd, 0o700)
    return fd

many = compile_mem("fm171_many", "ginibre_q3/character_ring_iter/fm39/sec166_class_pattern_flips.cpp")
flip = compile_mem("fm171_flip", "ginibre_q3/character_ring_iter/fm39/sec166_flip_single.cpp")
gp = compile_mem("fm171_gp", "ginibre_q3/character_ring_iter/fm39/sec162_gp_single.cpp",
                 openmp=False)

profiles = set()
def add(d):
    c = tuple(sorted((int(n), int(m)) for n, m in d.items() if m > 0))
    W = sum(n*m for n, m in c)
    if 48 <= W <= 200 and len(c) <= 9:
        profiles.add(c)

base = [2, 1, 4, 3, 4]
for t in range(1, 5):
    d = {i+1: t*base[i] for i in range(5)}
    add(d)
    for i in range(5):
        for change in (-1, 1):
            q = d.copy()
            q[i+1] += change
            add(q)
    for n in (6, 7, 8):
        q = d.copy()
        q[n] = 1
        add(q)

for K in range(8, 10):
    add({1: 2, **{n: 1 for n in range(2, K+1)}})
    add({1: 2, **{n: 2 for n in range(2, K+1)}})

fixed = [
 {1:3,2:1,3:4,4:3,5:4,6:1},
 {1:4,2:1,3:4,4:3,5:4,7:1},
 {1:2,2:3,3:4,4:3,5:4,6:2},
 {1:2,2:1,3:4,4:3,5:4,6:1,7:1},
 {1:2,2:1,3:4,4:3,5:4,6:1,7:1,8:1},
 {1:2,2:2,3:3,4:3,5:3,6:2},
 {1:3,2:4,3:4,4:3,5:3,6:2},
 {1:4,2:4,3:4,4:4,5:4},
 {1:2,2:3,3:6,4:5,5:4,6:3},
 {1:2,2:5,3:5,4:5,5:5,6:5},
 {1:2,2:4,3:4,4:4,5:4,6:4,7:4},
 {1:2,2:6,3:6,4:6,5:6},
 {1:2,2:5,3:4,4:5,5:4,6:5,7:4},
 {1:2,2:8,3:8,4:6,5:6},
 {1:2,2:2,3:8,4:8,5:8,6:4},
 {1:3,2:2,3:5,4:5,6:5,8:5},
 {1:2,2:4,3:5,5:6,6:6,7:6,8:4},
]
for d in fixed:
    add(d)

rng = random.Random(17120261002)
for _ in range(30):
    labels = sorted(rng.sample(range(2, 9), rng.randint(3, 6)))
    d = {1: rng.randint(2, 5)}
    for n in labels:
        d[n] = rng.randint(1, 6)
    add(d)

runs = []
for c in sorted(profiles, key=lambda z: (sum(n*m for n,m in z), z)):
    W = sum(n*m for n,m in c)
    mx = max(n for n,m in c)
    dmin = max(8, mx)
    pmin = max(6, mx)
    if (pmin-W) % 2:
        pmin += 1
    pmax = W - 2*dmin
    if (pmax-W) % 2:
        pmax -= 1
    if pmin > pmax:
        continue
    ps = {pmin, pmax, pmin+2*((pmax-pmin)//4),
          pmin+2*((pmax-pmin)//2)}
    if c == tuple((i+1, base[i]) for i in range(5)):
        ps.add(8)
    for p in sorted(ps):
        if pmin <= p <= pmax and (W-p) % 2 == 0:
            runs.append((c, p, W))

assert len(profiles) == 97 and len(runs) == 292
print("profiles", len(profiles), "profile_p_runs", len(runs), flush=True)

records = []
bands = {
    "48-99": [0, 0, 0, 0],
    "100-149": [0, 0, 0, 0],
    "150-200": [0, 0, 0, 0],
}
for i, (c, p, W) in enumerate(runs, 1):
    N = sum(m for _,m in c) + 1
    print(time.strftime("%H:%M:%S"), "start", i, "/", len(runs),
          "W", W, "p", p, "classes", len(c), "factors", N, flush=True)
    z = subprocess.run([f"/proc/self/fd/{many}", str(p),
                        *[f"{n}:{m}" for n,m in c]],
                       pass_fds=(many,), env=env, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if z.returncode:
        raise RuntimeError(z.stderr)
    key = "48-99" if W < 100 else ("100-149" if W < 150 else "150-200")
    stat = bands[key]
    stat[0] += 1
    for line in z.stdout.splitlines():
        if line.startswith("NOFLIP pattern"):
            q = re.search(r"minus-classes=(.*?) sigma p=([+-]?\d+) "
                          r"TopPair\(([-+]?\d+),([-+]?\d+)\) "
                          r"(ok|FAIL) g=(\d+) child=(\d+)", line)
            if not q:
                raise ValueError(line)
            minus, sp, a, b, status, parent, child = q.groups()
            delta = (W-p)//2
            ones = next((m for n,m in c if n == 1), 0)
            records.append((W, N, p, delta, ones, minus.strip(), int(sp),
                            int(a), int(b), status, int(parent), int(child),
                            Fraction(abs(int(a))+abs(int(b)), delta)))
        elif line.startswith("patterns="):
            q = re.search(r"patterns=(\d+) noflip=(\d+) TopPair_fail=(\d+)", line)
            pats, nf, fail = map(int, q.groups())
            stat[1] += pats
            stat[2] += nf
            stat[3] += fail
    print(time.strftime("%H:%M:%S"), "done", i, "/", len(runs),
          "patterns", stat[1], "no-flip records", len(records), flush=True)

print("BAND COUNTS", bands)
for r in records:
    print("NOFLIP", r)
assert sum(x[1] for x in bands.values()) == 17944
assert sum(x[2] for x in bands.values()) == 6
assert sum(x[3] for x in bands.values()) == 0

expected = {
 (48,8,"2 4", 5, 5,4075371,160741),
 (48,8,"1 2 3 4 5",-5,-5,4075371,160741),
 (50,6,"2 4", 5, 5,10625415,424915),
 (50,6,"1 2 3 4 5",-5,-5,10625415,424915),
 (56,8,"2 4",-4, 8,27576641,961895),
 (56,8,"1 2 3 4 5",-4,8,27576641,961895),
}
got = {(r[0],r[2],r[5],r[7],r[8],r[10],r[11]) for r in records}
assert got == expected
assert [bands[k][0] for k in ("48-99","100-149","150-200")] == [142,93,57]
print("full class-pattern census: PASS")

cases = [
 (8,[1,1,-2]+[3]*4+[-4]*3+[5]*4,[5,5],4075371,160741),
 (8,[-1,-1,-2]+[-3]*4+[-4]*3+[-5]*4,[-5,-5],4075371,160741),
 (6,[1,1]+[-2]*2+[3]*4+[-4]*3+[5]*4,[5,5],10625415,424915),
 (6,[-1,-1]+[-2]*2+[-3]*4+[-4]*3+[-5]*4,[-5,-5],10625415,424915),
 (8,[1,1,-2]+[3]*4+[-4]*3+[5]*4+[8],[-4,8],27576641,961895),
 (8,[-1,-1,-2]+[-3]*4+[-4]*3+[-5]*4+[8],[-4,8],27576641,961895),
]
for p,B,R,parent,child in cases:
    sp = -1 if sum(x < 0 for x in B) % 2 else 1
    q = subprocess.run([f"/proc/self/fd/{flip}", str(sp*p),
                        *map(str,B)], pass_fds=(flip,), env=env,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert q.returncode == 0 and "NO flip descent among" in q.stdout
    q = subprocess.run([f"/proc/self/fd/{gp}", str(p), *map(str,B),
                        "--", *map(str,R)], pass_fds=(gp,), env=env,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert q.returncode == 0
    assert f"g_p(B) = {parent}" in q.stdout
    assert f"g_p(B-R) = {child}" in q.stdout and "monotone" in q.stdout
print("six no-flip and TopPair cross-checks: PASS")
