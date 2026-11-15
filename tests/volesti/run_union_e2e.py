#!/usr/bin/env python3
"""End-to-end UNION volume accuracy suite (test plan T7) for ttc's LRA engine.

Runs a ttc binary on every tests/volesti/union_smt2/*.smt2 file with
--seed 1..K and compares the final `s volume <v>` line with the exact
volume in the file header:

    ; exact-volume: <number>      (required)
    ; family: <name>              (optional, used for the per-family summary)
    ; dim: <n>                    (optional, informational)
    ; timeout: <sec>              (optional, per-file minimum timeout)

Pass criterion per formula (Theorem: Pr[est in (1 +- eps) vol] >= 1 - delta):
  * failures = #seeds with est outside [(1-eps) vol, (1+eps) vol]
    (crashes, timeouts and missing `s volume` lines count as failures);
    PASS needs failures <= c_max, where c_max is the 1 - alpha quantile of
    Binomial(K, delta), alpha = 0.01 one-sided: if the true failure rate
    were <= delta, a count above c_max would occur with probability <= 1%.
  * |mean(est)/vol - 1| <= eps/2 (the mean relative error over the
    successful seeds; the estimator is meant to be (nearly) unbiased, so a
    mean off by eps/2 is a systematic error, not noise).
Column `z` = mean.err / (sd/sqrt(K)) is informational: |z| > 3 flags a
systematic bias that is real (not seed noise) even when it is inside eps/2.

eps/delta are ttc's defaults (0.8 / 0.2) unless --epsilon/--delta are given,
in which case they are also passed to ttc.

Usage:
    python3 tests/volesti/run_union_e2e.py [--ttc ./build/ttc] [--seeds 20]
        [--epsilon E] [--delta D] [--filter REGEX] [-j J] [--timeout SEC]
        [--extra-args "--nogmp"] [--json out.json] [--nice N] [--quiet]

Exit status is non-zero iff some formula FAILs.
"""
import argparse
import glob
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SMT2_DIR = os.path.join(HERE, "union_smt2")
DEFAULT_EPS = 0.8    # ttc --epsilon default
DEFAULT_DELTA = 0.2  # ttc --delta default
ALPHA = 0.01

VOL_RE = re.compile(r"^s volume\s+(\S+)", re.M)


def parse_header(path):
    meta = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*;\s*([\w-]+)\s*:\s*(.*?)\s*$", line)
            if not m:
                if line.strip() and not line.lstrip().startswith(";"):
                    break
                continue
            key, val = m.group(1).lower(), m.group(2)
            if key == "exact-volume":
                meta["exact"] = float(val)
            elif key == "timeout":
                meta["timeout"] = float(val)
            elif key == "dim":
                meta["dim"] = int(val)
            elif key == "family":
                meta["family"] = val
    meta.setdefault("family", os.path.basename(path).split("_")[0])
    meta.setdefault("dim", 0)
    return meta


def binom_cdf(c, n, p):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(c + 1))


def max_failures(n, delta, alpha=ALPHA):
    """Smallest c with P[Bin(n, delta) > c] <= alpha."""
    for c in range(n + 1):
        if 1.0 - binom_cdf(c, n, delta) <= alpha:
            return c
    return n


def run_one(ttc, path, seed, timeout, extra, nice):
    cmd = ([ "nice", "-n", str(nice)] if nice is not None else []) + [ttc, path, "--seed", str(seed)] + extra
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=timeout)
        out, rc = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired:
        return {"vol": None, "rc": None, "time": time.time() - t0, "note": "TIMEOUT"}
    dt = time.time() - t0
    m = VOL_RE.findall(out)
    vol = None
    if m:
        try:
            vol = float(m[-1])
        except ValueError:
            vol = None
    note = ""
    if rc != 0:
        note = "rc=%d" % rc
    elif vol is None:
        note = "no 's volume'"
    return {"vol": vol, "rc": rc, "time": dt, "note": note}


def summarize(name, meta, runs, eps, delta, K):
    exact = meta["exact"]
    ratios = [r["vol"] / exact for r in runs if r["vol"] is not None and r["rc"] == 0
              and math.isfinite(r["vol"])]
    bad = [r for r in runs if r["vol"] is None or r["rc"] != 0 or not math.isfinite(r["vol"])]
    fails = len(bad) + sum(1 for q in ratios if abs(q - 1) > eps)
    mean = sum(ratios) / len(ratios) - 1 if ratios else float("nan")
    sd = (sum((q - 1 - mean) ** 2 for q in ratios) / (len(ratios) - 1)) ** 0.5 if len(ratios) > 1 else float("nan")
    cmax = max_failures(K, delta)
    # Informational: z-score of the mean ratio against 1 (bias in units of its standard error).
    z = mean / (sd / math.sqrt(len(ratios))) if len(ratios) > 1 and sd > 0 else float("nan")
    ok = fails <= cmax and ratios and abs(mean) <= eps / 2
    return {"name": name, "family": meta["family"], "dim": meta["dim"], "exact": exact,
            "mean_rel_err": mean, "sd_rel_err": sd, "bias_z": z, "fails": fails, "cmax": cmax, "K": K,
            "crashes": len(bad), "notes": sorted(set(r["note"] for r in bad if r["note"])),
            "min_ratio": min(ratios) if ratios else None, "max_ratio": max(ratios) if ratios else None,
            "mean_time": sum(r["time"] for r in runs) / len(runs),
            "status": "PASS" if ok else "FAIL", "ests": [r["vol"] for r in runs]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ttc", default=os.path.join(REPO, "build", "ttc"))
    ap.add_argument("--seeds", type=int, default=20, help="run each file with --seed 1..K (default 20)")
    ap.add_argument("--epsilon", type=float, default=None, help="passed to ttc; default: ttc's 0.8")
    ap.add_argument("--delta", type=float, default=None, help="passed to ttc; default: ttc's 0.2")
    ap.add_argument("--filter", default=None, help="only files whose name (without .smt2) matches REGEX")
    ap.add_argument("--timeout", type=float, default=60.0, help="per-run timeout (s); header may raise it")
    ap.add_argument("--extra-args", default="", help="extra arguments passed to ttc, e.g. '--nogmp'")
    ap.add_argument("--jobs", "-j", type=int, default=2, help="parallel ttc runs (default 2)")
    ap.add_argument("--nice", type=int, default=19, help="run ttc under nice -n N (default 19; -1 disables)")
    ap.add_argument("--json", default=None, help="write per-formula results to this JSON file")
    ap.add_argument("--dir", default=SMT2_DIR, help="directory of .smt2 files")
    ap.add_argument("--quiet", action="store_true", help="only print the per-family summary")
    # Allow `--extra-args --nogmp` (argparse would read --nogmp as an option).
    argv = sys.argv[1:]
    for i, a in enumerate(argv[:-1]):
        if a == "--extra-args":
            argv[i:i + 2] = ["--extra-args=" + argv[i + 1]]
            break
    args = ap.parse_args(argv)

    ttc = os.path.abspath(args.ttc)
    if not os.access(ttc, os.X_OK):
        print("ttc binary not found/executable: %s" % ttc, file=sys.stderr)
        return 2
    files = sorted(glob.glob(os.path.join(args.dir, "*.smt2")))
    if args.filter:
        rx = re.compile(args.filter)
        files = [f for f in files if rx.search(os.path.basename(f)[:-5])]
    if not files:
        print("no test files selected", file=sys.stderr)
        return 2
    eps = DEFAULT_EPS if args.epsilon is None else args.epsilon
    delta = DEFAULT_DELTA if args.delta is None else args.delta
    extra = shlex.split(args.extra_args)
    if args.epsilon is not None:
        extra += ["--epsilon", repr(args.epsilon)]
    if args.delta is not None:
        extra += ["--delta", repr(args.delta)]
    K = args.seeds
    nice = None if args.nice < 0 else args.nice
    cmax = max_failures(K, delta)
    print("ttc: %s %s" % (ttc, " ".join(extra)))
    print("eps=%g delta=%g K=%d: PASS iff fails <= %d (Bin(K,delta) 1-%g quantile) and |mean rel err| <= %g"
          % (eps, delta, K, cmax, ALPHA, eps / 2))

    jobs = []
    metas = {}
    for path in files:
        meta = parse_header(path)
        if "exact" not in meta:
            print("%-28s (no exact-volume header, skipped)" % os.path.basename(path))
            continue
        metas[path] = meta
        for seed in range(1, K + 1):
            jobs.append((path, seed))
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = {(p, s): pool.submit(run_one, ttc, p, s, max(args.timeout, metas[p].get("timeout", 0)), extra, nice)
                for p, s in jobs}
        results = []
        hdr = "%-26s %-7s %3s %11s %9s %8s %6s %7s %7s %7s %7s  %s" % (
            "file", "family", "dim", "exact", "mean.err", "sd", "z", "fails", "minq", "maxq", "time", "status")
        if not args.quiet:
            print(hdr)
            print("-" * len(hdr))
        for path in metas:
            runs = [futs[(path, s)].result() for s in range(1, K + 1)]
            name = os.path.basename(path)[:-5]
            r = summarize(name, metas[path], runs, eps, delta, K)
            results.append(r)
            if not args.quiet:
                print("%-26s %-7s %3d %11.5g %+9.4f %8.4f %+6.1f %3d/%-3d %7s %7s %6.2fs  %s%s" % (
                    name[:26], r["family"][:7], r["dim"], r["exact"], r["mean_rel_err"], r["sd_rel_err"], r["bias_z"],
                    r["fails"], r["cmax"],
                    "-" if r["min_ratio"] is None else "%.3f" % r["min_ratio"],
                    "-" if r["max_ratio"] is None else "%.3f" % r["max_ratio"],
                    r["mean_time"], r["status"],
                    ("  [" + "; ".join(r["notes"]) + "]") if r["notes"] else ""), flush=True)

    # per-family summary
    fams = {}
    for r in results:
        fams.setdefault(r["family"], []).append(r)
    print()
    print("%-8s %5s %10s %10s %9s %9s %8s" % ("family", "files", "pass", "|mean.err|", "fails", "runs", "time"))
    for fam, rs in sorted(fams.items()):
        mae = [abs(r["mean_rel_err"]) for r in rs if r["mean_rel_err"] == r["mean_rel_err"]]
        print("%-8s %5d %10s %10.4f %9d %9d %7.2fs" % (
            fam, len(rs), "%d/%d" % (sum(r["status"] == "PASS" for r in rs), len(rs)),
            sum(mae) / len(mae) if mae else float("nan"), sum(r["fails"] for r in rs),
            sum(r["K"] for r in rs), sum(r["mean_time"] for r in rs) / len(rs)))
    n_fail = sum(r["status"] == "FAIL" for r in results)
    print("formulas: pass %d, fail %d; wall %.1fs" % (len(results) - n_fail, n_fail, time.time() - t_start))
    if n_fail:
        print("failed: " + ", ".join(r["name"] for r in results if r["status"] == "FAIL"))
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"ttc": ttc, "extra": extra, "eps": eps, "delta": delta, "K": K, "results": results}, f, indent=1)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
