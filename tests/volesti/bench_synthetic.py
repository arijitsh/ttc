#!/usr/bin/env python3
"""Benchmark two ttc binaries on random syntheticLRA instances.

For each instance both binaries run concurrently (same machine conditions),
niced, with a timeout.  Writes results.csv next to this script.
"""
import csv, lzma, os, random, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
INST = os.path.expanduser("~/code/codex/new-ttc-results/instances/syntheticLRA1131")
ARCHIVE = os.path.expanduser("~/code/codex/new-ttc-results/comparison/instances.csv")
BINS = {
    "start": sys.argv[1],   # update @ 75e05e0
    "now": sys.argv[2],     # current HEAD
}
N, SEED, TIMEOUT = 10, int(sys.argv[3]) if len(sys.argv) > 3 else 20261114, 900

files = sorted(f for f in os.listdir(INST) if f.endswith(".smt2.xz"))
pick = random.Random(SEED).sample(files, N)
# optional slice "i:j" of the picked instances (to run several benches in parallel)
SLICE = sys.argv[4] if len(sys.argv) > 4 else ":"
lo, hi = (int(v) if v else None for v in SLICE.split(":"))
pick = pick[lo:hi]
OUT = "results.csv" if SLICE == ":" else f"results_{SLICE.replace(':', '_')}.csv"
archive = {r["instance"]: r for r in csv.DictReader(open(ARCHIVE))}

def run(binary, path):
    t0 = time.perf_counter()
    try:
        p = subprocess.run(["nice", "-n", "19", binary, path, "-v", "2"], capture_output=True,
                           text=True, timeout=TIMEOUT)
        out, status = p.stdout + p.stderr, ("ok" if p.returncode == 0 else f"rc={p.returncode}")
    except subprocess.TimeoutExpired as e:
        out, status = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout"
    wall = time.perf_counter() - t0
    m = re.search(r"^s volume (\S+)", out, re.M)
    vol = float(m.group(1)) if m else float("nan")
    oracle = sum(float(x) for x in re.findall(r"polytope volume \S+ .* in (\S+)s", out))
    return wall, vol, status, oracle

rows = []
os.makedirs(os.path.join(HERE, "inst"), exist_ok=True)
for f in pick:
    name = f[:-3]
    path = os.path.join(HERE, "inst", name)
    if not os.path.exists(path):
        open(path, "w").write(lzma.open(os.path.join(INST, f), "rt").read())
    m = re.search(r"dim(\d+)_numpol(\d+)", name)
    with ThreadPoolExecutor(2) as ex:
        res = dict(zip(BINS, ex.map(lambda b: run(BINS[b], path), BINS)))
    a = archive.get(name, {})
    row = {"instance": name, "dim": m.group(1), "polytopes": m.group(2)}
    for k, (wall, vol, status, oracle) in res.items():
        row.update({f"{k}_wall": f"{wall:.2f}", f"{k}_volume": f"{vol:.6g}", f"{k}_status": status,
                    f"{k}_oracle_s": f"{oracle:.2f}"})
    row["archived_ttc_correct_volume"] = a.get("ttc_correct_volume", "")
    rows.append(row)
    print(row, flush=True)
    with open(os.path.join(HERE, OUT), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
