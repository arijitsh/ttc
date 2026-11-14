#!/usr/bin/env python3
"""Overnight job queue for ttc fuzzing / test campaigns.

night_runner.py --root DIR
  DIR/jobs/*.json   job specs {"name", "cmd" (shell string), "timeout" (s), "cwd"?}
                    picked up in name order; new files may be dropped in while
                    it runs (e.g. jobs added later in the night).
  DIR/workers.txt   optional: current worker limit (re-read before each launch)
  DIR/results/<name>.log    stdout+stderr of each job
  DIR/results/jobs.jsonl    one status line per finished job
  DIR/runner.log            runner progress
Worker limit: workers.txt if present, else 6 before 00:30 and 10 after
(the user's cap is 10).  Stops launching at the --stop time (default 08:00)
and exits when nothing is running and no job is left.
"""
import argparse
import datetime as dt
import json
import os
import signal
import subprocess
import time

CAP = 10


def now():
    return dt.datetime.now()


def worker_limit(root):
    try:
        with open(os.path.join(root, "workers.txt")) as f:
            return max(0, min(CAP, int(f.read().strip())))
    except (OSError, ValueError):
        t = now().time()
        return 10 if dt.time(0, 30) <= t < dt.time(12, 0) else 6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--stop", default="08:00")
    a = ap.parse_args()
    root = a.root
    jobs_dir = os.path.join(root, "jobs")
    res_dir = os.path.join(root, "results")
    os.makedirs(res_dir, exist_ok=True)
    hh, mm = map(int, a.stop.split(":"))
    stop_at = now().replace(hour=hh, minute=mm, second=0, microsecond=0)
    if stop_at <= now():
        stop_at += dt.timedelta(days=1)
    log = open(os.path.join(root, "runner.log"), "a", buffering=1)
    done = set()
    if os.path.exists(os.path.join(res_dir, "jobs.jsonl")):
        for line in open(os.path.join(res_dir, "jobs.jsonl")):
            done.add(json.loads(line)["name"])
    running = {}  # name -> (Popen, start, spec, logfile)
    log.write("%s runner start, stop launching at %s\n" % (now(), stop_at))
    while True:
        # reap
        for name, (p, t0, spec, lf) in list(running.items()):
            rc = p.poll()
            elapsed = time.time() - t0
            if rc is None and elapsed > spec.get("timeout", 4 * 3600):
                try:
                    os.killpg(p.pid, signal.SIGKILL)  # the shell and every child
                except ProcessLookupError:
                    pass
                p.wait()
                rc = "timeout"
            if rc is not None:
                lf.close()
                rec = {"name": name, "rc": rc, "start": t0, "wall": round(elapsed, 1),
                       "cmd": spec["cmd"]}
                with open(os.path.join(res_dir, "jobs.jsonl"), "a") as jf:
                    jf.write(json.dumps(rec) + "\n")
                log.write("%s done %s rc=%s %.0fs\n" % (now(), name, rc, elapsed))
                done.add(name)
                del running[name]
        pending = sorted(f[:-5] for f in os.listdir(jobs_dir)
                         if f.endswith(".json") and f[:-5] not in done and f[:-5] not in running)
        if now() < stop_at:
            while pending and len(running) < worker_limit(root):
                name = pending.pop(0)
                spec = json.load(open(os.path.join(jobs_dir, name + ".json")))
                lf = open(os.path.join(res_dir, name + ".log"), "w")
                p = subprocess.Popen(spec["cmd"], shell=True, stdout=lf, stderr=subprocess.STDOUT,
                                     cwd=spec.get("cwd", root), start_new_session=True)
                running[name] = (p, time.time(), spec, lf)
                log.write("%s start %s (running %d, limit %d)\n"
                          % (now(), name, len(running), worker_limit(root)))
        if not running and (not pending or now() >= stop_at):
            break
        time.sleep(5)
    log.write("%s runner exit (%d jobs done)\n" % (now(), len(done)))


if __name__ == "__main__":
    main()
