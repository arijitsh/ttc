#!/usr/bin/env python3
"""Summarize a night_runner.py campaign (see agents-todo-union-e2e.md).

night_summary.py --root NIGHTDIR [--eps 0.8]
Prints: job status counts; exact-volume union fuzz (bu) accuracy by variant
and dimension, failures and worst cases; syntheticLRA runs (time, volume vs
the archived ttc_correct volume); slow unit suites; e2e summary lines.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
from collections import defaultdict

ARCHIVE = os.path.expanduser("~/code/codex/new-ttc-results/comparison/instances.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--eps", type=float, default=0.8)
    a = ap.parse_args()
    R = a.root
    jobs = [json.loads(l) for l in open(os.path.join(R, "results", "jobs.jsonl"))] \
        if os.path.exists(os.path.join(R, "results", "jobs.jsonl")) else []
    st = defaultdict(int)
    for j in jobs:
        st[str(j["rc"])] += 1
    print("== jobs: %d finished, rc counts %s" % (len(jobs), dict(st)))
    for j in jobs:
        if j["rc"] not in (0,):
            print("   non-zero: %s rc=%s %.0fs" % (j["name"], j["rc"], j["wall"]))

    # exact-volume union fuzz
    cmd_of = {j["name"]: j["cmd"] for j in jobs}
    specs = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(os.path.join(R, "jobs", "*.json"))}
    variant_of_seed = {}
    for name, sp in specs.items():
        m = re.search(r"--seed (\d+).*--extra-args '([^']*)'", sp["cmd"])
        if "run_batch.py" in sp["cmd"] and m:
            variant_of_seed[int(m.group(1))] = m.group(2) or "default"
    rows = []
    for p in glob.glob(os.path.join(R, "bu", "batch_*.csv")):
        seed = int(re.search(r"batch_(\d+)\.csv", p).group(1))
        for r in csv.DictReader(open(p)):
            r["variant"] = variant_of_seed.get(seed, "?")
            rows.append(r)
    print("\n== exact-volume Boolean unions: %d runs" % len(rows))
    by = defaultdict(list)
    for r in rows:
        by[(r["variant"], int(r["dim"]))].append(r)
    print("%-28s %3s %5s %5s %6s %8s %8s %8s %6s" % ("variant", "dim", "runs", "ok", "fail",
                                                        "meanRel", "sdRel", "maxAbs", "tmed"))
    worst = []
    for (v, d), rs in sorted(by.items()):
        eps = 0.3 if "0.3" in v else a.eps
        ok = [r for r in rs if r["status"] == "ok"]
        rel = [float(r["rel_err"]) for r in ok]
        fail = sum(1 for x in rel if abs(x) > eps)
        mean = sum(rel) / len(rel) if rel else float("nan")
        sd = math.sqrt(sum((x - mean) ** 2 for x in rel) / max(1, len(rel) - 1)) if rel else float("nan")
        mx = max((abs(x) for x in rel), default=float("nan"))
        walls = sorted(float(r["wall"]) for r in rs)
        print("%-28s %3d %5d %5d %6d %+8.4f %8.4f %8.4f %6.1f" % (
            v[:28], d, len(rs), len(ok), fail, mean, sd, mx, walls[len(walls) // 2] if walls else 0))
        for r in rs:
            if r["status"] != "ok" or abs(float(r["rel_err"])) > eps:
                worst.append((v, r))
    print("\n   non-ok or |rel err| > eps (first 40):")
    for v, r in worst[:40]:
        print("   %-24s %s seed %s: %s vol=%s exact=%s rel=%s" % (v[:24], r["file"], r["seed"], r["status"],
                                                                  r["volume"], r["exact"], r["rel_err"]))

    # syntheticLRA
    arch = {}
    if os.path.exists(ARCHIVE):
        for r in csv.DictReader(open(ARCHIVE)):
            arch[r["instance"]] = r
    print("\n== syntheticLRA sample")
    print("%-52s %8s %12s %12s %7s %9s" % ("instance", "wall", "volume", "archived", "ratio", "arch_wall"))
    for name, sp in sorted(specs.items()):
        if "_syn_" not in name:
            continue
        log = os.path.join(R, "results", name + ".log")
        inst = re.search(r"/syn/(\S+\.smt2)", sp["cmd"]).group(1)
        txt = open(log).read() if os.path.exists(log) else ""
        m = re.search(r"^s volume (\S+)", txt, re.M)
        w = re.search(r"wall ([\d.]+) s", txt)
        vol = float(m.group(1)) if m else float("nan")
        ar = arch.get(inst, {})
        av = float(ar["ttc_correct_volume"]) if ar.get("ttc_correct_volume") else float("nan")
        print("%-52s %8s %12.5g %12.5g %7.3f %9s" % (inst[:52], w.group(1) if w else "-", vol, av,
                                                    vol / av if av == av and av else float("nan"),
                                                    ar.get("ttc_correct_wall_seconds", "-")))

    print("\n== slow unit suites / e2e")
    for name in sorted(specs):
        if "_slow_" in name or "_e2e_" in name or "fuzz" in name or "union_e2e" in name:
            log = os.path.join(R, "results", name + ".log")
            txt = open(log).read() if os.path.exists(log) else "(not run)"
            tail = [l for l in txt.splitlines() if re.search(r"test cases|Status|^pass |FAIL|TOTAL", l)]
            print("%-40s %s" % (name, " | ".join(tail[-3:]) if tail else txt.strip().splitlines()[-1:]))


if __name__ == "__main__":
    main()
