#!/usr/bin/env python3
import sys, os, subprocess, re

DUMP_DIR = "Datasets/Gwfa256"

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
        "run_args": [DUMP_DIR, "762"],
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

if __name__ == "__main__":
    main()
