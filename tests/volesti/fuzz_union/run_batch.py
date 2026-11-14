#!/usr/bin/env python3
"""Generate one batch of exact-volume union formulas and run ttc on them.

run_batch.py --ttc BIN --work DIR --seed S [--count N] [--seeds K] [--timeout T]
             [--dims 2,3,4] [--extra-args "..."] -- [generator args]
Writes DIR/batch_S.csv with one row per (formula, ttc seed):
  file,dim,exact,seed,status,volume,rel_err,wall
status: ok | timeout | error(rc) | noout
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ttc", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--count", type=int, default=25)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=180)
    ap.add_argument("--dims", default="2,3,4")
    ap.add_argument("--extra-args", default="")
    ap.add_argument("gen_args", nargs="*")
    a = ap.parse_args()
    fdir = os.path.join(a.work, "batch_%d" % a.seed)
    subprocess.run([sys.executable, os.path.join(HERE, "gen_bool_unions.py"), "--out", fdir,
                    "--count", str(a.count), "--seed", str(a.seed), "--dims", a.dims]
                   + a.gen_args, check=True)
    out = os.path.join(a.work, "batch_%d.csv" % a.seed)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "dim", "exact", "seed", "status", "volume", "rel_err", "wall"])
        for f in sorted(os.listdir(fdir)):
            path = os.path.join(fdir, f)
            head = open(path).read(400)
            exact = float(re.search(r"exact-volume: (\S+)", head).group(1))
            dim = int(re.search(r"; dim: (\d+)", head).group(1))
            for s in range(1, a.seeds + 1):
                cmd = ["nice", "-n", "19", a.ttc, path, "--seed", str(1000 * a.seed + s)]
                cmd += a.extra_args.split()
                t0 = time.time()
                try:
                    p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
                    m = re.search(r"^s volume (\S+)", p.stdout + p.stderr, re.M)
                    status = "ok" if (p.returncode == 0 and m) else (
                        "error(%d)" % p.returncode if p.returncode else "noout")
                    vol = float(m.group(1)) if m else float("nan")
                except subprocess.TimeoutExpired:
                    status, vol = "timeout", float("nan")
                wall = time.time() - t0
                rel = vol / exact - 1 if exact and vol == vol else float("nan")
                w.writerow([f, dim, exact, s, status, vol, rel, round(wall, 3)])
                fh.flush()
    print("wrote", out)


if __name__ == "__main__":
    main()
