#!/usr/bin/env python3
"""SMT-LIB syntax fuzzing of ttc's QF_LRA volume engine.

Generates (tests/volesti/fuzz_lra/gen_fuzz.py) tiny satisfiable, bounded QF_LRA
formulas that exercise SMT-LIB syntax, runs ttc on each one and classifies:

  OK           `s volume v` with v within tolerance of a Monte-Carlo estimate of
               the evaluator (|v/mc - 1| <= --rel-tol, or |v - mc| <= --z * se)
  WRONG        finite `s volume` far from the MC estimate
  CRASH        signal / abort / non-zero exit / no `s volume` and no `Error:`
  UNSUPPORTED  explicit `Error: ...` message
  TIMEOUT      no answer within --timeout seconds

Prints counts per class, per syntax feature, and per failure signature (class +
normalised message, with the features over-represented in that signature).
Failing formulas are shrunk (commands and subterms deleted while the failure
signature persists; a plain bounding-box assert is kept so the MC reference
stays valid) and saved under tests/volesti/fuzz_lra/failures/.

ttc runs under `nice -n 19`; at most --jobs (default 2) CPU-bound tasks (ttc
processes or MC evaluations) run at once.

Usage: python3 tests/volesti/run_fuzz_lra.py [--count 500] [--seed 1]
          [--ttc build-gcc/ttc] [--timeout 20] [--jobs 2] [--regen]
          [--no-shrink] [--shrink-per-sig 2] [--filter RE]
"""
import argparse
import collections
import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
import threading
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FUZZ = os.path.join(HERE, "fuzz_lra")
sys.path.insert(0, FUZZ)
sys.path.insert(0, HERE)

import gen_fuzz  # noqa: E402

# Features whose failures share a root cause found by this fuzzer (see the report).
MASK_DEFAULT = ("define0_real,define0_bool,define_args_real,define_args_bool,"
                "define_param_shadow,define_nested,real_ite,to_real_in_qflra")
CLASSES = ["OK", "WRONG", "CRASH", "UNSUPPORTED", "TIMEOUT"]
CPU = threading.Semaphore(2)

# ---------------------------------------------------------------------------
# ttc + classification
# ---------------------------------------------------------------------------


def run_ttc(ttc, path, timeout):
    t0 = time.time()
    with CPU:
        try:
            r = subprocess.run(["nice", "-n", "19", ttc, path], capture_output=True,
                               text=True, errors="replace", timeout=timeout)
            out, rc, to = r.stdout + r.stderr, r.returncode, False
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or b"").decode("utf-8", "replace") if isinstance(e.stdout, bytes) \
                else (e.stdout or "")
            rc, to = None, True
    return {"out": out, "rc": rc, "timeout": to, "time": time.time() - t0}


def normalise(msg):
    msg = re.sub(r"[-+]?\d+(\.\d+)?(e[-+]?\d+)?", "N", msg)
    return msg.strip()[:140]


def classify(run, ref, rel_tol, ztol):
    """(class, message, volume)."""
    out = run["out"]
    if run["timeout"]:
        m = re.search(r"^c polytopes: (\d+)", out, re.M)
        return "TIMEOUT", "timeout (%s)" % ("polytopes %s" % m.group(1) if m else
                                            "before the polytope count"), None
    m = re.search(r"^s volume (\S+)", out, re.M)
    vol = None
    if m:
        try:
            vol = float(m.group(1))
        except ValueError:
            vol = None
    if m and run["rc"] == 0 and vol is not None:
        mc, se = ref
        if vol != vol or vol in (float("inf"), float("-inf")):
            return "WRONG", "non-finite volume %s" % m.group(1), vol
        if abs(vol - mc) <= max(rel_tol * mc, ztol * se):
            return "OK", "", vol
        kind = "zero" if vol == 0 else ("low" if vol < mc else "high")
        return "WRONG", "volume %s vs MC (%s)" % (kind, "ratio<1" if vol < mc else "ratio>1"), vol
    what = re.search(r"what\(\):\s*(.*)", out)
    err = re.search(r"^(Error:.*|.*[Ee]rror:.*)$", out, re.M)
    rc = run["rc"]
    if rc is not None and (rc < 0 or rc >= 128):
        msg = what.group(1) if what else (err.group(1) if err else "")
        return "CRASH", "rc=%d %s" % (rc, msg), vol
    if err:
        return "UNSUPPORTED", err.group(1).strip(), vol
    if not m:
        sline = re.search(r"^s .*$", out, re.M)
        return "CRASH", "rc=%s no 's volume' (%s)" % (rc, sline.group(0) if sline else
                                                       "no s-line"), vol
    return "CRASH", "rc=%s" % rc, vol


def signature(cls, msg):
    return "%s: %s" % (cls, normalise(msg))


# ---------------------------------------------------------------------------
# S-expressions that keep the raw tokens (quoted symbols, strings)
# ---------------------------------------------------------------------------
TOKEN_RE = re.compile(r'\s+|;[^\n]*|\(|\)|\|[^|]*\||"(?:[^"]|"")*"|[^\s();|"]+')


def parse_raw(text):
    stack = [[]]
    for m in TOKEN_RE.finditer(text):
        tok = m.group(0)
        if tok[0].isspace() or tok[0] == ";":
            continue
        if tok == "(":
            stack.append([])
        elif tok == ")":
            done = stack.pop()
            stack[-1].append(done)
        else:
            stack[-1].append(tok)
    return stack[0]


def render(e):
    if isinstance(e, list):
        return "(" + " ".join(render(x) for x in e) + ")"
    return e


def render_cmds(cmds, meta_line):
    return meta_line + "\n" + "\n".join(render(c) for c in cmds) + "\n"


def qconst(q):
    q = Fraction(q)
    a = abs(q)
    s = str(a.numerator) if a.denominator == 1 else "(/ %d %d)" % (a.numerator, a.denominator)
    return "(- %s)" % s if q < 0 else s


def canonical_bbox(meta):
    parts = []
    for v, (lo, hi) in zip(meta["vars"], meta["box"]):
        parts.append("(<= %s %s)" % (qconst(lo), v))
        parts.append("(<= %s %s)" % (v, qconst(hi)))
    return parse_raw("(assert (and %s))" % " ".join(parts))[0]


def insert_bbox(cmds, box_cmd):
    heads = [c[0] if isinstance(c, list) and c else None for c in cmds]
    pos = len(cmds)
    for i in range(len(cmds) - 1, -1, -1):
        if heads[i] == "check-sat":
            pos = i
            break
    else:
        if "exit" in heads:
            pos = heads.index("exit")
    return cmds[:pos] + [box_cmd] + cmds[pos:]


def get_at(e, path):
    for k in path:
        e = e[k]
    return e


def set_at(e, path, val):
    if not path:
        return val
    out = list(e)
    out[path[0]] = set_at(e[path[0]], path[1:], val)
    return out


def node_paths(e, path=()):
    if isinstance(e, list):
        yield path, e
        for i, x in enumerate(e):
            yield from node_paths(x, path + (i,))


NARY = {"and", "or", "+", "-", "*", "=", "distinct", "<", "<=", ">", ">=", "=>", "xor"}
BOOLISH = {"and", "or", "not", "=>", "xor", "<", "<=", ">", ">=", "=", "distinct", "!", "let",
           "ite"}


def candidates(cmds, protected):
    """Smaller variants of the command list, biggest deletions first."""
    for i in range(len(cmds) - 1, -1, -1):
        if i in protected:
            continue
        yield cmds[:i] + cmds[i + 1:]
    for ci, c in enumerate(cmds):
        if ci in protected or not isinstance(c, list) or not c:
            continue
        if c[0] not in ("assert", "define-fun"):
            continue
        for path, node in node_paths(c):
            if not path or not node or isinstance(node[0], list):
                continue
            head = node[0]
            if head in ("let",) and len(node) == 3:
                yield cmds[:ci] + [set_at(c, path, node[2])] + cmds[ci + 1:]
                if len(node[1]) > 1:
                    for j in range(len(node[1])):
                        b = node[1][:j] + node[1][j + 1:]
                        yield cmds[:ci] + [set_at(c, path, [head, b, node[2]])] + cmds[ci + 1:]
                continue
            if head == "!":
                yield cmds[:ci] + [set_at(c, path, node[1])] + cmds[ci + 1:]
                continue
            if head in ("and", "or", "not", "ite", "=>", "xor"):
                for j in range(1, len(node)):
                    yield cmds[:ci] + [set_at(c, path, node[j])] + cmds[ci + 1:]
            if head in NARY and len(node) > 3:
                for j in range(1, len(node)):
                    yield cmds[:ci] + [set_at(c, path, node[:j] + node[j + 1:])] + cmds[ci + 1:]
            if head in BOOLISH:
                yield cmds[:ci] + [set_at(c, path, "true")] + cmds[ci + 1:]


def shrink(text, meta, target_sig, check, budget):
    """Greedy delta-debugging on commands and subterms."""
    meta_line = text.splitlines()[0]
    cmds = parse_raw(text)
    box = canonical_bbox(meta)
    with_box = insert_bbox(cmds, box)
    used = [0]

    def ok(cs):
        if used[0] >= budget:
            return False
        used[0] += 1
        return check(render_cmds(cs, meta_line)) == target_sig

    if ok(with_box):
        cmds = with_box
    elif not ok(cmds):
        return text, used[0]                        # re-rendering lost the failure
    while used[0] < budget:
        prot = {i for i, c in enumerate(cmds) if c is box or
            (isinstance(c, list) and c and c[0] in ("declare-const", "declare-fun"))}
        for cand in candidates(cmds, prot):
            if ok(cand):
                cmds = cand
                break
        else:
            break
    return render_cmds(cmds, meta_line), used[0]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ttc", default=os.path.join(REPO, "build-gcc", "ttc"))
    ap.add_argument("--count", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--gen-dir", default=os.path.join(FUZZ, "gen"))
    ap.add_argument("--fail-dir", default=os.path.join(FUZZ, "failures"))
    ap.add_argument("--json", default=os.path.join(FUZZ, "last_run.json"))
    ap.add_argument("--regen", action="store_true", help="regenerate the corpus")
    ap.add_argument("--timeout", type=float, default=20.0)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--mc-points", type=int, default=200000)
    ap.add_argument("--rel-tol", type=float, default=0.25)
    ap.add_argument("--z", type=float, default=4.0)
    ap.add_argument("--filter", default=None, help="regex on file names")
    ap.add_argument("--mask", default=MASK_DEFAULT,
                    help="features with an already-known root cause; a second per-feature "
                         "table covers only the formulas using none of them")
    ap.add_argument("--no-shrink", action="store_true")
    ap.add_argument("--shrink-per-sig", type=int, default=2)
    ap.add_argument("--shrink-budget", type=int, default=250, help="ttc runs per shrink")
    args = ap.parse_args()
    global CPU
    CPU = threading.Semaphore(args.jobs)
    try:
        os.nice(19)
    except OSError:
        pass

    # ---- corpus ------------------------------------------------------------
    os.makedirs(args.gen_dir, exist_ok=True)
    paths = []
    for idx in range(args.count):
        p = os.path.join(args.gen_dir, "fuzz_%04d.smt2" % idx)
        if args.regen or not os.path.exists(p):
            text, _ = gen_fuzz.generate_one(idx, args.seed, 0.03)
            with open(p, "w") as f:
                f.write(text)
        paths.append(p)
    if args.filter:
        paths = [p for p in paths if re.search(args.filter, os.path.basename(p))]
    items = []
    for p in paths:
        text = open(p).read()
        meta = gen_fuzz.read_meta(text)
        vol, se, frac = gen_fuzz.mc_volume(text, meta["box"], args.mc_points,
                                           seed=10007 + meta["index"])
        items.append({"path": p, "text": text, "meta": meta, "mc": vol, "se": se,
                      "frac": frac})
    print("corpus: %d formulas, MC references with %d points each" % (len(items),
                                                                      args.mc_points))
    sys.stdout.flush()

    # ---- run ttc -----------------------------------------------------------
    def job(it):
        run = run_ttc(args.ttc, it["path"], args.timeout)
        cls, msg, vol = classify(run, (it["mc"], it["se"]), args.rel_tol, args.z)
        it.update(cls=cls, msg=msg, vol=vol, time=run["time"],
                  tail=run["out"].strip().splitlines()[-6:], sig=signature(cls, msg))
        return it

    t0 = time.time()
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        for k, it in enumerate(ex.map(job, items)):
            if (k + 1) % 50 == 0:
                print("  %d/%d done (%.0fs)" % (k + 1, len(items), time.time() - t0))
                sys.stdout.flush()

    # ---- report ------------------------------------------------------------
    cnt = collections.Counter(it["cls"] for it in items)
    print("\nclasses: " + "  ".join("%s=%d" % (c, cnt[c]) for c in CLASSES))
    feats = sorted({f for it in items for f in it["meta"]["features"]})
    mask = set(filter(None, args.mask.split(",")))
    clean = [it for it in items if not mask & set(it["meta"]["features"])]
    cc = collections.Counter(it["cls"] for it in clean)
    print("without masked features (%s): %d formulas, " % (",".join(sorted(mask)) or "-",
                                                           len(clean))
          + "  ".join("%s=%d" % (c, cc[c]) for c in CLASSES))
    hdr = "%-22s %5s %5s %5s %5s %5s %5s" % ("feature", "n", *[c[:5] for c in CLASSES])
    print("\n%s   | unmasked: %s" % (hdr, hdr[23:]))
    for f in feats:
        row = []
        for pool in (items, clean):
            sub = [it for it in pool if f in it["meta"]["features"]]
            c = collections.Counter(it["cls"] for it in sub)
            row.append("%5d %5d %5d %5d %5d %5d" % (len(sub), *[c[k] for k in CLASSES]))
        print("%-22s %s   |           %s" % (f, row[0], row[1]))
    sigs = collections.defaultdict(list)
    for it in items:
        if it["cls"] != "OK":
            sigs[it["sig"]].append(it)
    base = {f: sum(f in it["meta"]["features"] for it in items) / len(items) for f in feats}
    print("\nfailure signatures:")
    for s, its in sorted(sigs.items(), key=lambda kv: -len(kv[1])):
        fr = {f: sum(f in it["meta"]["features"] for it in its) / len(its) for f in feats}
        always = [f for f in feats if fr[f] == 1.0 and base[f] < 0.9]
        lift = sorted((f for f in feats if fr[f] >= 0.5 and base[f] < 0.9),
                      key=lambda f: -fr[f] / max(base[f], 1e-9))[:4]
        print("  %4d  %s" % (len(its), s))
        print("        in all: %s" % (", ".join(always) or "-"))
        print("        enriched: %s" % ", ".join("%s(%.0f%%/%.0f%%)" % (f, 100 * fr[f],
                                                                     100 * base[f]) for f in lift))
        print("        e.g. %s" % ", ".join(os.path.basename(it["path"]) for it in its[:5]))

    # ---- shrink --------------------------------------------------------------
    if not args.no_shrink and sigs:
        os.makedirs(args.fail_dir, exist_ok=True)
        work = os.path.join(args.fail_dir, ".work")
        os.makedirs(work, exist_ok=True)
        jobs = []
        for s, its in sigs.items():
            its = sorted(its, key=lambda it: len(it["text"]))
            jobs += [(s, it) for it in its[:args.shrink_per_sig]]

        def do_shrink(job_):
            s, it = job_
            meta = it["meta"]
            tmp = os.path.join(work, "cand_%d_%d.smt2" % (meta["index"], threading.get_ident()))
            cls = it["cls"]
            budget = 0 if cls == "TIMEOUT" else args.shrink_budget

            def check(text):
                with open(tmp, "w") as f:
                    f.write(text)
                ref = (it["mc"], it["se"])
                if cls == "WRONG":
                    try:
                        with CPU:
                            vol, se, frac = gen_fuzz.mc_volume(text, meta["box"], 50000,
                                                               seed=meta["index"])
                    except Exception:
                        return None
                    if frac < 0.01:
                        return None
                    ref = (vol, se)
                run = run_ttc(args.ttc, tmp, min(args.timeout, 10))
                c, m, _ = classify(run, ref, args.rel_tol, args.z)
                return signature(c, m)

            small, used = shrink(it["text"], meta, s, check, budget) if budget else (it["text"], 0)
            slug = re.sub(r"[^A-Za-z0-9]+", "_", s)[:60].strip("_")
            out = os.path.join(args.fail_dir, "%s__%04d.smt2" % (slug, meta["index"]))
            with open(out, "w") as f:
                f.write("; class: %s\n; signature: %s\n; original: %s (%d bytes, %d shrink runs)\n"
                        % (cls, s, os.path.relpath(it["path"], REPO), len(it["text"]), used))
                f.write("; MC volume (original): %.6g +- %.2g, ttc: %s\n" % (
                    it["mc"], it["se"], it["vol"]))
                for line in it["tail"][-4:]:
                    f.write(";   ttc> %s\n" % line[:200])
                f.write(small)
            os.remove(tmp) if os.path.exists(tmp) else None
            return s, out, len(small)

        print("\nshrinking %d failing formulas ..." % len(jobs))
        sys.stdout.flush()
        with cf.ThreadPoolExecutor(args.jobs) as ex:
            for s, out, size in ex.map(do_shrink, jobs):
                print("  %-70s -> %s (%d bytes)" % (s[:70], os.path.relpath(out, REPO), size))
                sys.stdout.flush()
        try:
            os.rmdir(work)
        except OSError:
            pass

    if args.json:
        with open(args.json, "w") as f:
            json.dump([{k: it[k] for k in ("path", "cls", "msg", "vol", "mc", "se", "time",
                                           "tail")} | {"features": it["meta"]["features"]}
                       for it in items], f, indent=1)
    print("\nclasses: " + "  ".join("%s=%d" % (c, cnt[c]) for c in CLASSES))
    return 0 if cnt["OK"] == len(items) else 1


if __name__ == "__main__":
    sys.exit(main())
