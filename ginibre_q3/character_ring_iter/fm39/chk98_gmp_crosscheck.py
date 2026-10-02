import os, re, subprocess

ROOT = "/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39/"
ENV = dict(os.environ, TMPDIR="/dev/shm", OMP_NUM_THREADS="16")

def compile_memfd(source_name):
    with open(ROOT + source_name, "r") as f:
        source = f.read()
    fd = os.memfd_create("fmchk98_" + source_name.replace(".", "_"), flags=0)
    cmd = ["g++", "-pipe", "-std=c++17", "-O2", "-fopenmp",
           "-x", "c++", "-o", f"/proc/self/fd/{fd}", "-",
           "-lgmpxx", "-lgmp"]
    cp = subprocess.run(cmd, input=source, text=True, capture_output=True,
                        pass_fds=(fd,), env=ENV)
    assert cp.returncode == 0, cp.stderr
    os.set_inheritable(fd, True)
    return fd

def run(fd, args):
    return subprocess.run([f"/proc/self/fd/{fd}", *args], text=True,
                          capture_output=True, pass_fds=(fd,), env=ENV)

fd = compile_memfd("sec166_flip_single.cpp")
out = run(fd, ["62", "-40", "42", "-44", "46", "48",
               "50", "52", "54", "56", "58", "60"])
os.close(fd)
assert out.returncode == 0, out.stderr
signs = re.findall(r"pair \d+/66 .* sign ([+-])", out.stdout)
assert len(signs) == 66 and all(s == "-" for s in signs)
assert "NO flip descent among 66 pairs" in out.stdout
print("GMP sec166_flip_single: 66/66 negative; no flip descent PASS")

fd = compile_memfd("sec162_gp_single.cpp")
out = run(fd, ["62", "-40", "42", "-44", "46", "48",
               "50", "52", "54", "56", "58", "60", "--", "58", "60"])
os.close(fd)
assert out.returncode == 0, out.stderr
assert "g_p(B) = 453207534222864" in out.stdout
assert "g_p(B-R) = 179646349605  monotone" in out.stdout
print("GMP sec162_gp_single: parent/child values PASS")

fd = compile_memfd("sec162_gp_single.cpp")
args = ["6"] + ["-1"]*13 + ["-2"] + ["-3"]*5 + ["-4"] + ["--", "-2", "-4"]
out = run(fd, args)
os.close(fd)
assert out.returncode == 0, out.stderr
assert "g_p(B) = 700607800" in out.stdout
assert "g_p(B-R) = 723462700" in out.stdout
print("GMP sec162_gp_single: ratio-3/7 parent/child values PASS")

fd = compile_memfd("sec166_class_pattern_flips.cpp")
out = run(fd, ["6", "1:13", "2:1", "3:5", "4:1"])
os.close(fd)
assert out.returncode == 0, out.stderr
assert "profile W=34 p=6 delta=14 classes=4 factors=21 residual=yes" in out.stdout
assert "patterns=16 noflip=0 TopPair_fail=0" in out.stdout
print("GMP sec166_class_pattern_flips: 16 patterns, all have flip descent PASS")
