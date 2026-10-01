import os, re, sys, time, random, subprocess, collections

if any(x in sys.argv[1:] for x in ("-h", "--help")):
    print("Run the FM-SEC167 seeded search and consecutive-label checks.")
    raise SystemExit(0)

ROOT = "/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39/"

def build(src, name, flags, libs):
    fd = os.memfd_create(name, 0)
    os.set_inheritable(fd, True)
    z = subprocess.run(
        ["g++", "-std=c++17", "-O3", "-pipe"] + flags
        + ["-x", "c++", "-", "-o", f"/proc/self/fd/{fd}"] + libs,
        input=src.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        pass_fds=(fd,), env={**os.environ, "TMPDIR": "/dev/shm"})
    if z.returncode:
        raise RuntimeError(z.stderr.decode()[-3000:])
    return fd

src = open(ROOT + "sec166_flip_single.cpp").read()
src = src.replace(
    "vector<int> Lam=B;Lam.push_back(p);",
    "bool hill = getenv(\"SEC167_HILL\") != nullptr;\n"
    "  vector<int> Lam=B;Lam.push_back(p);")
old = """    int a=vals[prs[k].first],b=vals[prs[k].second];vector<int> C=Lam;C.erase(find(C.begin(),C.end(),a));C.erase(find(C.begin(),C.end(),b));"""
new = """    int a=vals[prs[k].first],b=vals[prs[k].second];
    if(hill){
      bool ap=(abs(a)==abs(p)),bp=(abs(b)==abs(p));
      if(ap&&bp)continue;
      if(ap||bp){
        int x=ap?b:a;
        if(count(B.begin(),B.end(),x)!=1)continue;
      }else if(a!=b&&
               (count(B.begin(),B.end(),a)!=1 ||
                count(B.begin(),B.end(),b)!=1))continue;
    }
    vector<int> C=Lam;
    C.erase(find(C.begin(),C.end(),a));
    C.erase(find(C.begin(),C.end(),b));"""
assert old in src
src = src.replace(old, new)
src = src.replace(
    'if(!found)printf("NO flip descent among %zu pairs\\\\n",prs.size());',
    'if(!found){if(hill)printf("NO pair-free B-lowering flip\\\\n");'
    'else printf("NO flip descent among %zu pairs\\\\n",prs.size());}')

flipfd = build(src, "sec167_flip", ["-fopenmp"], ["-lgmpxx", "-lgmp"])
gpfd = build(open(ROOT + "sec162_gp_single.cpp").read(),
             "sec167_gp", [], ["-lgmpxx", "-lgmp"])
FLIP = f"/proc/self/fd/{flipfd}"
GP = f"/proc/self/fd/{gpfd}"
ENV = {**os.environ, "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"}

def flip(p, B, hill):
    env = ENV.copy()
    if hill:
        env["SEC167_HILL"] = "1"
    else:
        env.pop("SEC167_HILL", None)
    z = subprocess.run(
        [FLIP, str(p)] + list(map(str, B)),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env=env, pass_fds=(flipfd,), timeout=60)
    if z.returncode:
        raise RuntimeError(z.stderr + "\n" + z.stdout[-500:])
    return z.stdout

def gp(p, B, removals):
    args = [GP, str(p)] + list(map(str, B))
    for R in removals:
        args += ["--"] + list(map(str, R))
    z = subprocess.run(
        args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env=ENV, pass_fds=(gpfd,), timeout=60)
    if z.returncode:
        raise RuntimeError(z.stderr + "\n" + z.stdout[-500:])
    vals = []
    for line in z.stdout.splitlines():
        a = re.search(r"g_p\(B\) = (-?\d+)", line)
        b = re.search(r"g_p\(B-R\) = (-?\d+)", line)
        if a:
            vals.append(int(a.group(1)))
        elif b:
            vals.append(int(b.group(1)))
    return z.stdout, vals

def make_profile(rng, lo, hi, typ):
    for _ in range(800):
        if typ == 0:
            B = [rng.randint(1, 9) if rng.random() < .78
                 else rng.randint(10, 24)
                 for _ in range(rng.randint(12, 34))]
        elif typ == 1:
            k = rng.randint(6, 20)
            step = 1 if rng.random() < .65 else 2
            B = []
            for a in range(1 if step == 1 else 3, k + 1, step):
                B += [a] * rng.randint(1, 3)
        elif typ == 2:
            k = rng.randint(10, 28)
            B = []
            for a in range(1, k + 1):
                B += [a] * max(0, k // a - 1)
            B = B[:55]
        else:
            B = []
            for a in range(1, rng.randint(8, 17)):
                B += [a] * rng.randint(1, 6)

        W = sum(B)
        if not lo <= W <= hi or sum(a >= 3 for a in B) < 2:
            continue
        signs = {a: (1 if rng.randrange(2) else -1) for a in set(B)}
        B = [a * signs[a] for a in B]
        mx = max(map(abs, B))
        ps = [p for p in range(max(6, mx + 1), min(W - 16, W - 2 * mx) + 1)
              if (W - p) % 2 == 0]
        if ps:
            return B, rng.choice(ps), W
    return None

def apply_flip(B, a, b, p):
    C = B.copy()
    if abs(a) == p or abs(b) == p:
        x = b if abs(a) == p else a
        C.remove(x)
        C.append(-x)
    elif a == b:
        C.remove(a)
        C.remove(a)
        C += [-a, -a]
    else:
        C.remove(a)
        C.remove(b)
        C += [-a, -b]
    return C

def top_pair(B):
    choices = []
    for i, a in enumerate(B):
        for b in B[i + 1:]:
            if (abs(a) - abs(b)) % 2 == 0:
                choices.append(((abs(a) + abs(b), max(abs(a), abs(b)),
                                 min(abs(a), abs(b))), (a, b)))
    return max(choices, key=lambda z: (z[0], tuple(sorted(z[1]))))[1] \
        if choices else None

def all_even_removals(B):
    ct = collections.Counter(B)
    out = set()
    for x in ct:
        if abs(x) % 2 == 0:
            out.add((x,))
    vals = list(ct)
    for i, a in enumerate(vals):
        for b in vals[i:]:
            if (abs(a) - abs(b)) % 2 == 0 and (a != b or ct[a] >= 2):
                out.add(tuple(sorted((a, b))))
    return sorted(out)

rng = random.Random(1672031)
bands = [(41, 80), (81, 120), (121, 200), (201, 300)]
stats = collections.defaultdict(collections.Counter)
for band in bands:
    for q in range(6):
        item = make_profile(rng, *band, q % 4)
        if item is None:
            stats[band]["generation_reject"] += 1
            continue
        B, p, W = item
        steps = 0
        status = "cap"
        for steps in range(100):
            sigma = -1 if sum(x < 0 for x in B) % 2 else 1
            out = flip(sigma * p, B, True)
            m = re.search(
                r"FLIP DESCENT at \(([+-]?\d+),([+-]?\d+)\): "
                r"eps_v\*A = (-?\d+)", out)
            if m:
                a, b, value = map(int, m.groups())
                if value == 0:
                    status = "zero_preserving"
                    break
                B = apply_flip(B, a, b, p)
                continue
            if "NO pair-free B-lowering flip" in out:
                full = flip(sigma * p, B, False)
                status = ("true_noflip" if "NO flip descent among" in full
                          else "other_nonnegative_flip")
                break
        stats[band]["trials"] += 1
        stats[band][status] += 1
        stats[band]["steps"] += steps
        print("PROFILE", band, "W", W, "steps", steps,
              "terminal", status, flush=True)

        if status == "true_noflip":
            tp = top_pair(B)
            even = max((x for x in B if abs(x) % 2 == 0),
                       key=lambda x: abs(x), default=None)
            Rs = [tp] + ([(even,)] if even is not None and (even,) != tp else [])
            _, values = gp(p, B, Rs)
            parent = values[0]
            assert parent >= values[1]
            if even is not None:
                assert parent >= values[2]
            if not (parent >= values[1] or
                    (even is not None and parent >= values[2])):
                _, all_values = gp(p, B, all_even_removals(B))
                assert any(parent >= x for x in all_values[1:])

print("SUMMARY")
for band in bands:
    print(band, dict(stats[band]))
print("TOTAL", sum(v["trials"] for v in stats.values()),
      "true_noflip", sum(v["true_noflip"] for v in stats.values()))

expected = {
    11: 97, 12: 17182, 15: 1126500, 16: 75512019,
    19: 35048990879, 20: 2622754205406,
    23: 3617099112665551, 24: 235839418360029549
}
for k in range(8, 25):
    B = [-n if n % 4 == 1 else n for n in range(1, k)]
    W = k * (k - 1) // 2
    if (W - k) % 2:
        continue
    delta = (W - k) // 2
    if delta < 8 or max(map(abs, B)) > delta:
        continue
    sigma = -1 if sum(x < 0 for x in B) % 2 else 1
    out = flip(sigma * k, B, False)
    m = re.search(
        r"FLIP DESCENT at \(([+-]?\d+),([+-]?\d+)\): "
        r"eps_v\*A = (-?\d+)", out)
    if m:
        a, b, value = map(int, m.groups())
        assert value == expected[k]
        print("CONSEC", k, W, delta, "flip", a, b, value)
    else:
        assert "NO flip descent among" in out and k == 8
        tp = top_pair(B)
        even = max(x for x in B if abs(x) % 2 == 0)
        _, values = gp(k, B, [tp, (even,)])
        assert values == [478, 18, 88]
        print("CONSEC", k, W, delta, "no-flip", values,
              "TopPair", tp, "TopEven", even)

# Independent polynomial/Catalan evaluation of the k=8 control.
def catalan(n):
    from math import comb
    return comb(2*n, n) // (n + 1)

def character(n):
    if n == 0:
        return [1]
    if n == 1:
        return [0, 1]
    a, b = [1], [0, 1]
    for _ in range(1, n):
        c = [0] * (len(b) + 1)
        for i, v in enumerate(b):
            c[i + 1] += v
        for i, v in enumerate(a):
            c[i] -= v
        a, b = b, c
    return b

def direct_g(p, B):
    F = {(0, 0): 1}
    for z in B:
        u = character(abs(z))
        eps = 1 if z > 0 else -1
        factor = collections.defaultdict(int)
        for i, coeff in enumerate(u):
            if coeff:
                factor[i, 0] += coeff
                factor[0, i] += eps * coeff
        nxt = collections.defaultdict(int)
        for (i, j), x in F.items():
            for (k, ell), y in factor.items():
                nxt[i + k, j + ell] += x * y
        F = nxt
    up = character(p)
    ans = 0
    for (i, j), x in F.items():
        if j % 2:
            continue
        for k, coeff in enumerate(up):
            if coeff and (i + k) % 2 == 0:
                ans += x * coeff * catalan(j // 2) * catalan((i + k) // 2)
    return ans

B = [-1, 2, 3, 4, -5, 6, 7]
assert [direct_g(8, B),
        direct_g(8, [x for x in B if x not in (-5, 7)]),
        direct_g(8, [x for x in B if x != 6])] == [478, 18, 88]
print("CATALAN_CONTROL", 478, 18, 88)
print("PASS")
