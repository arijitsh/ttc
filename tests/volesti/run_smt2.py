#!/usr/bin/env python3
"""End-to-end single-polytope volume suite for ttc's LRA volume engine.

Runs ttc on every tests/volesti/smt2/*.smt2 file, parses the final
`s volume <v>` line and compares it with the exact volume given in the
file header:

    ; exact-volume: <number>
    ; max-rel-error: <tol>      (optional, default 0.2)
    ; dim: <n>                  (optional, informational)
    ; report-only: true         (optional: print the result, never FAIL)
    ; timeout: <sec>            (optional: per-file minimum timeout)

Usage:
    python3 tests/volesti/run_smt2.py [--ttc ./build/ttc] [--seeds K]
                                      [--filter REGEX] [--timeout SEC]
                                      [--extra-args "..."] [--jobs J]

Columns: `poly` is the raw single-polytope oracle estimate (`polytope
volume` line = volume_cooling_gaussians output BEFORE ttc clamps it to
[(1-eps')(2r/sqrt n)^n, (1+eps')(2R)^n]); `got` is the final `s volume`
answer (clamped oracle value fed through the union estimator, which for a
single polytope multiplies it by ~1.00-1.01); rel.err is for `got`.  A note
`clamped x<ratio>` flags runs where got/poly is far from 1, i.e. ttc's
clamp (not the oracle) produced the answer.

With --seeds K each file is run with --seed 1..K; a file passes only if
every seed passes.  Exit status is non-zero iff any non-report-only case
fails (wrong volume, missing `s volume` line, crash, or timeout).
"""
import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SMT2_DIR = os.path.join(HERE, "smt2")
DEFAULT_TOL = 0.2

# Files whose header we do not own / that predate the header convention.
FALLBACK_META = {
    "long_box_2x1000.smt2": {"exact": 2000.0, "dim": 2},
}

VOL_RE = re.compile(r"^s volume\s+(\S+)", re.M)
POLY_RE = re.compile(r"^polytope volume\s+(\S+)", re.M)
WARN_RE = re.compile(r"(warn|error|exception|abort|terminate|nan|inf\b)", re.I)


def parse_header(path):
    meta = dict(FALLBACK_META.get(os.path.basename(path), {}))
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*;\s*([\w-]+)\s*:\s*(.*?)\s*$", line)
            if not m:
                continue
            key, val = m.group(1).lower(), m.group(2)
            if key == "exact-volume":
                meta["exact"] = float(val)
            elif key == "max-rel-error":
                meta["tol"] = float(val)
            elif key == "timeout":
                meta["timeout"] = float(val)
            elif key == "dim":
                meta["dim"] = int(val)
            elif key == "report-only":
                meta["report_only"] = val.lower() in ("1", "true", "yes")
    if "dim" not in meta:
        with open(path) as f:
            meta["dim"] = len(re.findall(r"\(declare-(?:const|fun)\s+\S+\s+(?:\(\s*\)\s*)?Real\)", f.read()))
    meta.setdefault("tol", DEFAULT_TOL)
    meta.setdefault("report_only", False)
    return meta


def run_one(ttc, path, seed, timeout, extra):
    cmd = [ttc, path, "-v", "2"] + extra
    if seed is not None:
        cmd += ["--seed", str(seed)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=timeout)
        out, rc = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        return {"vol": None, "poly": None, "rc": None, "time": time.time() - t0,
                "note": "TIMEOUT", "out": out}
    dt = time.time() - t0
    m = VOL_RE.findall(out)
    vol = float(m[-1]) if m else None
    poly = [float(x) for x in POLY_RE.findall(out)]
    notes = []
    if rc != 0:
        notes.append("rc=%d" % rc)
    if vol is None:
        notes.append("no 's volume'")
    # Interesting non-comment lines (warnings/errors), excluding known-benign ones.
    for line in out.splitlines():
        if line.startswith("c ") or line.startswith("s volume") or line.startswith("polytope volume"):
            continue
        if WARN_RE.search(line):
            notes.append(line.strip()[:80])
    return {"vol": vol, "poly": poly, "rc": rc, "time": dt, "note": "; ".join(notes), "out": out}


def rel_err(got, exact):
    if got is None:
        return None
    if exact == 0:
        return abs(got)
    return abs(got - exact) / abs(exact)


def fmt(x):
    if x is None:
        return "-"
    return "%.5g" % x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ttc", default=os.path.join(REPO, "build", "ttc"))
    ap.add_argument("--seeds", type=int, default=0,
                    help="run each file with --seed 1..K (default: ttc's default seed only)")
    ap.add_argument("--filter", default=None, help="only run files whose basename matches REGEX")
    ap.add_argument("--timeout", type=float, default=30.0, help="per-run timeout in seconds")
    ap.add_argument("--extra-args", default="", help="extra arguments passed to ttc")
    ap.add_argument("--jobs", "-j", type=int, default=min(8, os.cpu_count() or 1),
                    help="number of ttc runs in parallel (timings get noisier with more)")
    ap.add_argument("--show-output", action="store_true", help="dump ttc output for failing runs")
    args = ap.parse_args()

    ttc = os.path.abspath(args.ttc)
    if not os.access(ttc, os.X_OK):
        print("ttc binary not found/executable: %s" % ttc, file=sys.stderr)
        return 2
    files = sorted(glob.glob(os.path.join(SMT2_DIR, "*.smt2")))
    if args.filter:
        rx = re.compile(args.filter)
        files = [f for f in files if rx.search(os.path.basename(f))]
    if not files:
        print("no test files selected", file=sys.stderr)
        return 2
    seeds = list(range(1, args.seeds + 1)) if args.seeds > 0 else [None]
    extra = shlex.split(args.extra_args)

    hdr = "%-28s %3s %4s %11s %11s %11s %9s %7s  %-6s %s" % (
        "file", "dim", "seed", "exact", "poly", "got", "rel.err", "time", "status", "notes")
    print(hdr)
    print("-" * len(hdr))
    n_fail = n_pass = n_report = 0
    failed = []
    total_t0 = time.time()
    jobs = []
    for path in files:
        meta = parse_header(path)
        if "exact" not in meta:
            print("%-28s  (no exact-volume header, skipped)" % os.path.basename(path))
            continue
        for seed in seeds:
            jobs.append((path, meta, seed))
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = [pool.submit(run_one, ttc, path, seed, max(args.timeout, meta.get("timeout", 0)), extra) for path, meta, seed in jobs]
        for (path, meta, seed), fut in zip(jobs, futs):
            name = os.path.basename(path)
            r = fut.result()
            err = rel_err(r["vol"], meta["exact"])
            crashed = r["rc"] is None or r["rc"] != 0 or r["vol"] is None
            if meta["report_only"]:
                status = "REPORT"
                n_report += 1
            elif crashed or err is None or err > meta["tol"] or err != err:
                status = "FAIL"
                n_fail += 1
                failed.append((name, seed))
            else:
                status = "PASS"
                n_pass += 1
            poly = r["poly"][-1] if r["poly"] else None
            if poly and r["vol"] and poly > 0 and abs(r["vol"] / poly - 1) > 0.05:
                r["note"] = ("clamped x%.3g; " % (r["vol"] / poly)) + r["note"]
            print("%-28s %3d %4s %11s %11s %11s %9s %6.2fs  %-6s %s" % (
                name[:-5][:28], meta["dim"], "-" if seed is None else seed, fmt(meta["exact"]),
                fmt(poly), fmt(r["vol"]), fmt(err), r["time"], status, r["note"]), flush=True)
            if status == "FAIL" and args.show_output:
                print("    | " + "\n    | ".join(r["out"].splitlines()[-40:]))
    print("-" * len(hdr))
    print("pass %d, fail %d, report-only %d, wall %.1fs" % (n_pass, n_fail, n_report, time.time() - total_t0))
    if failed:
        print("failed: " + ", ".join("%s%s" % (n, "" if s is None else "@seed%d" % s) for n, s in failed))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
