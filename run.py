#!/usr/bin/env python3
import sys, os, subprocess, re

DUMP_DIR = "Datasets/Gwfa295"
GOLDEN_SCORES = os.path.join(
    DUMP_DIR, "trueScores.txt")
GOLDEN_DEBUG = os.path.join(
    DUMP_DIR, "wfDebugTrue0.txt")

TESTS = {
    "1": {
        "name": "fast",
        "make": ["make", "-j"],
        "run_args": [DUMP_DIR, "15"],
    },
    "2": {
        "name": "full",
        "make": ["make", "-j"],
        "run_args": [DUMP_DIR],
    },
    "3": {
        "name": "wfDebug",
        "make": ["make", "DBG=1", "-j"],
        "run_args": [DUMP_DIR, "1"],
    },
}

def run_cmd(cmd, capture=False):
    if capture:
        lines = []
        p = subprocess.Popen(cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            print(line, end="", flush=True)
            lines.append(line)
        p.wait()
        return p.returncode, "".join(lines)
    r = subprocess.run(cmd)
    return r.returncode, ""

def extract_kernel_time(output):
    m = re.search(r"kernel time:\s*(\d+)us",
        output)
    return int(m.group(1)) if m else None

def check_scores(scores_file, golden_file, limit=None):
    with open(scores_file) as f:
        sim = [l.strip() for l in f if l.strip()]
    with open(golden_file) as f:
        gold = [l.strip() for l in f if l.strip()]
    if limit and limit > 0:
        sim = sim[:limit]
        gold = gold[:limit]
    passed = failed = 0
    for i, (s, g) in enumerate(zip(sim, gold)):
        if s == g:
            passed += 1
        else:
            failed += 1
            print(f"  [{i}] FAIL sim={s} gold={g}")
    if len(sim) != len(gold):
        print(f"  count mismatch: "
              f"sim={len(sim)} gold={len(gold)}")
        failed += 1
    print(f"Scores: {passed} passed, "
          f"{failed} failed")
    return failed == 0

def check_debug(debug_file, golden_file):
    with open(debug_file) as f:
        sim = f.read()
    with open(golden_file) as f:
        gold = f.read()
    if sim == gold:
        print("Debug trace: PASS")
        return True
    else:
        print("Debug trace: FAIL (diff)")
        return False

def main():
    if len(sys.argv) != 2 \
            or sys.argv[1] not in TESTS:
        print("Usage: python3 run.py <1|2|3>")
        print("  1: fast (15 iterations)")
        print("  2: full (all iterations)")
        print("  3: wfDebug")
        sys.exit(1)

    t = TESTS[sys.argv[1]]
    print(f"=== {t['name']} ===")

    # clean + build
    for f in ["scores.txt", "wfDebug.txt"]:
        if os.path.exists(f):
            os.remove(f)
    run_cmd(["make", "clean"])
    rc, _ = run_cmd(t["make"])
    if rc != 0:
        print("Build failed"); sys.exit(1)

    # run kernel
    cmd = ["./gwfa"] + t["run_args"]
    rc, output = run_cmd(cmd, capture=True)
    if rc != 0:
        print(f"Kernel exited with code {rc}")
        sys.exit(1)

    # report kernel time
    kt = extract_kernel_time(output)
    if kt is not None:
        print(f"kernel time: {kt}us")

    # correctness check
    mode = sys.argv[1]
    ok = True
    if mode in ("1", "2"):
        limit = 15 if mode == "1" else None
        ok = check_scores(
            "scores.txt", GOLDEN_SCORES, limit)
    elif mode == "3":
        ok = check_scores(
            "scores.txt", GOLDEN_SCORES, 1)
        ok2 = check_debug(
            "wfDebug.txt", GOLDEN_DEBUG)
        ok = ok and ok2

    if not ok:
        sys.exit(1)

if __name__ == "__main__":
    main()
