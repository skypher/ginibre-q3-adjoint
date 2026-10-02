from collections import defaultdict
import ast
import gzip
import os
import pathlib
import random
import re
import shlex
import subprocess

def fuse(a, b):
    return range(abs(a-b), a+b+1, 2)

def product(F, G):
    H = defaultdict(int)
    for (a, b), x in F.items():
        for (c, d), y in G.items():
            for i in fuse(a, c):
                for j in fuse(b, d):
                    H[i, j] += x*y
    return {k: v for k, v in H.items() if v}

d = {(1, 0): 1, (0, 1): -1}
d2 = product(d, d)
chi11 = {(0, 0): 1, (1, 1): 1}
d2chi = product(d2, chi11)
prefix = sum(x*d2chi.get(key, 0) for key, x in d2.items()
             if sum(key) <= 2)
assert d2 == {(0, 0): 2, (2, 0): 1, (1, 1): -2, (0, 2): 1}
assert prefix == -6
print("PASS exact unrestricted-pair prefix T=2:", prefix, flush=True)

D = pathlib.Path(
    "/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39"
)
tree = ast.parse((D / "str7e_interior_budget.py").read_text())
cpp = next(
    ast.literal_eval(node.value)
    for node in ast.walk(tree)
    if isinstance(node, ast.Assign)
    and any(isinstance(x, ast.Name) and x.id == "CPP"
            for x in node.targets)
)
cpp = cpp.replace(
    "assert(jobs.size()==5531);",
    "assert(jobs.size()==1523);"
)
assert "assert(jobs.size()==1523);" in cpp

pat = re.compile(
    r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+) D="
)
jobs = []
for line in gzip.open(D / "sec166_census_w40_noflip.log.gz", "rt"):
    match = pat.match(line)
    if not match:
        continue
    W, p, B, phi = match.groups()
    word = tuple(map(int, B.split())) + (int(p),)
    if sum(x < 0 for x in word) != 4:
        continue
    assert sum(abs(x) for x in word[:-1]) == int(W)
    jobs.append((0, word, int(phi)))
assert len(jobs) == 1139

rng = random.Random(861593)
seen = {word for _, word, _ in jobs}

def candidate(core, plus):
    B = tuple(-x for x in core) + tuple(plus)
    W = sum(abs(x) for x in B)
    mx = max(map(abs, B))
    p = max(6, mx + 1, max(plus, default=0))
    if p % 2 != W % 2:
        p += 1
    if W > 80 or (W-p)//2 < max(8, mx):
        return None
    assert sum(x < 0 for x in B) == 4
    assert len(set(core)) == 4 and not (set(core) & set(plus))
    assert p > max(core)
    return B + (p,)

random_jobs = []
while len(random_jobs) < 256:
    core = tuple(sorted(rng.sample(range(3, 21), 4)))
    choices = [x for x in range(5, 21) if x not in core]
    plus = tuple(rng.choices(choices, k=rng.randint(4, 5)))
    word = candidate(core, plus)
    if word is not None and word not in seen:
        seen.add(word)
        random_jobs.append((1, word, -1))
jobs.extend(random_jobs)

fixed_jobs = []
core = (1, 2, 3, 4)
while len(fixed_jobs) < 128:
    plus = tuple(rng.choices(range(5, 21), k=rng.randint(4, 6)))
    word = candidate(core, plus)
    if word is not None and word not in seen:
        seen.add(word)
        fixed_jobs.append((2, word, -1))
jobs.extend(fixed_jobs)

assert len(jobs) == 1523
assert all(sum(x < 0 for x in word) == 4 for _, word, _ in jobs)
data = "".join(
    f"{group} {len(word)} {phi} " + " ".join(map(str, word)) + "\n"
    for group, word, phi in jobs
)
print(
    f"PREPARED census=1139 random_general={len(random_jobs)} "
    f"random_core1234={len(fixed_jobs)} total={len(jobs)}",
    flush=True
)

os.environ["TMPDIR"] = "/dev/shm"
obj = os.memfd_create("fmstr8f_obj", 0)
exe = os.memfd_create("fmstr8f_exe", 0)
subprocess.run(
    ["g++", "-Werror=return-type", "-O2", "-std=c++17", "-fopenmp",
     "-pipe", "-x", "c++", "-", "-c", "-o", f"/proc/self/fd/{obj}"],
    input=cpp, text=True, pass_fds=(obj,), check=True
)
link_info = subprocess.run(
    ["g++", "-###", "-fno-use-linker-plugin", "-fopenmp",
     f"/proc/self/fd/{obj}", "-o", f"/proc/self/fd/{exe}"],
    capture_output=True, text=True, check=True
)
link = next(
    shlex.split(line) for line in link_info.stderr.splitlines()
    if "/collect2 " in line
)
link[0] = "/usr/bin/ld"
subprocess.run(link, pass_fds=(obj, exe), check=True)
subprocess.run(
    [f"/proc/self/fd/{exe}", "--threads", "24"],
    input=data, text=True, pass_fds=(exe,), check=True
)
