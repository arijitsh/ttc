#!/usr/bin/env python3
"""Compare TTC weighted --PB against weighted_to_unweighted.py + ApproxMC.

The script consumes MCComp Track4 PWMC weighted CNF instances, converts each
instance to SMT-LIB with declare-projvar/declare-weight commands, runs TTC, then
runs the Python weighted-to-unweighted converter and ApproxMC on the original
CNF.  The ApproxMC result is scaled by the converter multiplier/divisor and
compared against TTC's weighted count.
"""

from __future__ import annotations

import argparse
import decimal
import gzip
import lzma
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Iterable


DEFAULT_REPO_URL = "https://github.com/arijitsh/mccomp-test-instances.git"
DEFAULT_TRACK_DIR = "Track4_PWMC"
DEFAULT_W2U = Path("/home/arijit/solvers/weighted-to-unweighted/weighted_to_unweighted.py")
APPROXMC_FALLBACKS = [
    Path("/home/arijit/solvers/sharpSMT/bin/approxmc"),
    Path("/home/arijit/solvers/sharpSMT/bins/approxmc"),
]


decimal.getcontext().prec = 80
Decimal = decimal.Decimal


@dataclass
class WeightedCnf:
    num_vars: int
    num_clauses: int
    clauses: list[list[int]]
    sampling_vars: list[int]
    weights: dict[int, Decimal]


@dataclass
class RunResult:
    instance: Path
    ttc_wmc: Decimal | None
    approxmc_wmc: Decimal | None
    rel_error: Decimal | None
    status: str
    detail: str


@dataclass
class PreparedReference:
    path: Path | None
    unsat: bool
    multiplier: Decimal


def run_cmd(cmd: list[str], timeout: int, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )


def resolve_executable(cli_value: str | None, fallbacks: Iterable[Path]) -> Path | str | None:
    if cli_value:
        return Path(cli_value)
    for candidate in fallbacks:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
    return None


def ensure_instances(args: argparse.Namespace, work_dir: Path) -> Path:
    if args.instances:
        return Path(args.instances).resolve()

    checkout = work_dir / "mccomp-test-instances"
    track = checkout / args.track_dir
    if track.exists():
        return track

    if not checkout.exists():
        cmd = [
            "git",
            "clone",
            "--depth",
            "1",
            "--filter=blob:none",
            "--sparse",
            args.repo_url,
            str(checkout),
        ]
        proc = run_cmd(cmd, args.timeout)
        if proc.returncode != 0:
            raise RuntimeError(
                "failed to clone instances repo\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
    proc = run_cmd(["git", "-C", str(checkout), "sparse-checkout", "set", args.track_dir], args.timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            "failed to configure sparse checkout\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    if not track.exists():
        raise RuntimeError(f"track directory not found after clone: {track}")
    return track


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if path.suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_decimal(token: str) -> Decimal:
    token = token.strip()
    if "/" in token:
        num, den = token.split("/", 1)
        return Decimal(num) / Decimal(den)
    return Decimal(token)


def parse_weighted_cnf(path: Path) -> WeightedCnf:
    num_vars: int | None = None
    num_clauses: int | None = None
    clauses: list[list[int]] = []
    sampling_vars: list[int] = []
    sampling_seen: set[int] = set()
    weights: dict[int, Decimal] = {}
    found_sampling_set = False

    with open_text(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) >= 4 and fields[0] == "p" and fields[1] == "cnf":
                num_vars = int(fields[2])
                num_clauses = int(fields[3])
                continue
            if fields[0] == "c":
                if len(fields) >= 4 and fields[1] == "p" and fields[2] == "show":
                    found_sampling_set = True
                    for tok in fields[3:]:
                        var = int(tok)
                        if var == 0:
                            break
                        if var > 0 and var not in sampling_seen:
                            sampling_seen.add(var)
                            sampling_vars.append(var)
                    continue
                if len(fields) >= 3 and fields[1] == "ind":
                    found_sampling_set = True
                    for tok in fields[2:]:
                        var = int(tok)
                        if var == 0:
                            break
                        if var > 0 and var not in sampling_seen:
                            sampling_seen.add(var)
                            sampling_vars.append(var)
                    continue
                if len(fields) >= 5 and fields[1] == "p" and fields[2] == "weight":
                    lit = int(fields[3])
                    if lit != 0:
                        weights[lit] = parse_decimal(fields[4])
                    continue
                continue
            if fields[0] == "w" and len(fields) >= 3:
                lit = int(fields[1])
                if lit != 0:
                    weights[lit] = parse_decimal(fields[2])
                continue
            if fields[0].lstrip("-").isdigit():
                clause: list[int] = []
                for tok in fields:
                    lit = int(tok)
                    if lit == 0:
                        break
                    clause.append(lit)
                clauses.append(clause)

    if num_vars is None or num_clauses is None:
        raise RuntimeError(f"missing p cnf header in {path}")
    if not found_sampling_set:
        sampling_vars = list(range(1, num_vars + 1))
    return WeightedCnf(num_vars, num_clauses, clauses, sampling_vars, weights)


def smt_name(var: int) -> str:
    return f"x{var}"


def literal_to_smt(lit: int) -> str:
    name = smt_name(abs(lit))
    return name if lit > 0 else f"(not {name})"


def cnf_to_smt2(cnf: WeightedCnf, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        out.write("(set-logic QF_UF)\n")
        for var in range(1, cnf.num_vars + 1):
            out.write(f"(declare-const {smt_name(var)} Bool)\n")

        chunk_size = 256
        for start in range(0, len(cnf.sampling_vars), chunk_size):
            chunk = cnf.sampling_vars[start : start + chunk_size]
            names = " ".join(smt_name(v) for v in chunk)
            out.write(f"(declare-projvar {names})\n")

        for lit, weight in sorted(cnf.weights.items(), key=lambda item: (abs(item[0]), item[0] < 0)):
            var = smt_name(abs(lit))
            if lit < 0:
                out.write(f"(declare-weight -{var} {weight})\n")
            else:
                out.write(f"(declare-weight {var} {weight})\n")

        for clause in cnf.clauses:
            if not clause:
                out.write("(assert false)\n")
            elif len(clause) == 1:
                out.write(f"(assert {literal_to_smt(clause[0])})\n")
            else:
                lits = " ".join(literal_to_smt(lit) for lit in clause)
                out.write(f"(assert (or {lits}))\n")
        out.write("(check-sat)\n")


def parse_ttc_wmc(stdout: str) -> Decimal:
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "s" and fields[1] in {"wmc", "mc"}:
            return parse_decimal(fields[2])
    raise RuntimeError(f"could not parse TTC count from output:\n{stdout}")


def parse_w2u_divisor(stdout: str) -> int:
    m = re.search(r"divide by:\s*2\*\*(\d+)", stdout)
    if not m:
        raise RuntimeError(f"could not parse weighted-to-unweighted divisor:\n{stdout}")
    return int(m.group(1))


def parse_w2u_multiplier(cnf_path: Path) -> Decimal:
    with cnf_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if "c MUST MULTIPLY BY" in line:
                fields = line.split()
                return parse_decimal(fields[4])
    return Decimal(1)


def normalize_projection_comments_for_approxmc(cnf_path: Path) -> None:
    text = cnf_path.read_text(encoding="utf-8", errors="replace")
    normalized = re.sub(r"^c p show\b", "c ind", text, flags=re.MULTILINE)
    if normalized != text:
        cnf_path.write_text(normalized, encoding="utf-8")


def parse_approxmc_count(stdout: str) -> Decimal:
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "s" and fields[1] == "mc":
            return Decimal(fields[2])
    m = re.search(r"Number of solutions is:\s*([0-9]+)\s*\*\s*2\*\*(\d+)", stdout)
    if m:
        return Decimal(m.group(1)) * (Decimal(2) ** int(m.group(2)))
    m = re.search(r"Number of solutions is:\s*([0-9]+)", stdout)
    if m:
        return Decimal(m.group(1))
    raise RuntimeError(f"could not parse ApproxMC count from output:\n{stdout}")


def relative_error(a: Decimal, b: Decimal) -> Decimal:
    denom = max(abs(a), abs(b), Decimal(1))
    return abs(a - b) / denom


def safe_artifact_name(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".xz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def declared_literal_weight(cnf: WeightedCnf, lit: int) -> Decimal:
    if lit in cnf.weights:
        return cnf.weights[lit]
    if -lit in cnf.weights:
        return Decimal(1) - cnf.weights[-lit]
    return Decimal(1)


def prepare_reference_cnf(cnf: WeightedCnf, work_dir: Path, artifact: str) -> PreparedReference:
    clauses = [list(clause) for clause in cnf.clauses]
    assignments: dict[int, bool] = {}
    queue = [clause[0] for clause in clauses if len(clause) == 1]

    while queue:
        lit = queue.pop()
        var = abs(lit)
        value = lit > 0
        previous = assignments.get(var)
        if previous is not None:
            if previous != value:
                return PreparedReference(None, True, Decimal(0))
            continue

        assignments[var] = value
        simplified_clauses: list[list[int]] = []
        for clause in clauses:
            simplified: list[int] = []
            satisfied = False
            for clause_lit in clause:
                if abs(clause_lit) != var:
                    simplified.append(clause_lit)
                    continue
                if (clause_lit > 0) == value:
                    satisfied = True
                    break
            if satisfied:
                continue
            if not simplified:
                return PreparedReference(None, True, Decimal(0))
            if len(simplified) == 1:
                queue.append(simplified[0])
            simplified_clauses.append(simplified)
        clauses = simplified_clauses

    sampling_set = set(cnf.sampling_vars)
    multiplier = Decimal(1)
    for var, value in assignments.items():
        if var in sampling_set:
            multiplier *= declared_literal_weight(cnf, var if value else -var)

    remaining_sampling = [var for var in cnf.sampling_vars if var not in assignments]
    remaining_weights = {
        lit: weight for lit, weight in cnf.weights.items() if abs(lit) not in assignments
    }

    output = work_dir / "reference" / f"{artifact}.preprocessed.cnf"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        out.write(f"p cnf {cnf.num_vars} {len(clauses)}\n")
        out.write("c p show")
        for var in remaining_sampling:
            out.write(f" {var}")
        out.write(" 0\n")
        for lit, weight in sorted(remaining_weights.items(), key=lambda item: (abs(item[0]), item[0] < 0)):
            out.write(f"w {lit} {weight} 0\n")
        for clause in clauses:
            out.write(" ".join(str(lit) for lit in clause))
            out.write(" 0\n")
        out.write(f"c MUST MULTIPLY BY {multiplier} 0\n")

    return PreparedReference(output, False, multiplier)


def run_one(instance: Path, args: argparse.Namespace, work_dir: Path) -> RunResult:
    artifact = safe_artifact_name(instance)
    smt_path = work_dir / "smt2" / f"{artifact}.smt2"
    converted_cnf = work_dir / "w2u" / f"{artifact}.unweighted.cnf"
    converted_cnf.parent.mkdir(parents=True, exist_ok=True)

    try:
        cnf = parse_weighted_cnf(instance)
        cnf_to_smt2(cnf, smt_path)
        reference_cnf = prepare_reference_cnf(cnf, work_dir, artifact)

        ttc_cmd = [
            str(args.ttc),
            "--PB",
            "-v",
            "0",
            "--seed",
            str(args.seed),
            str(smt_path),
        ]
        ttc_proc = run_cmd(ttc_cmd, args.timeout)
        if ttc_proc.returncode != 0:
            return RunResult(instance, None, None, None, "TTC_FAIL", ttc_proc.stderr.strip())
        ttc_wmc = parse_ttc_wmc(ttc_proc.stdout)

        if reference_cnf.unsat or reference_cnf.multiplier == 0:
            approx_wmc = Decimal(0)
            rel = relative_error(ttc_wmc, approx_wmc)
            status = "OK" if rel <= args.tolerance else "MISMATCH"
            detail = f"smt={smt_path} reference=unit-propagated-unsat-or-zero-weight"
            return RunResult(instance, ttc_wmc, approx_wmc, rel, status, detail)

        w2u_cmd = [
            sys.executable,
            str(args.w2u),
            "--prec",
            str(args.precision),
            str(reference_cnf.path),
            str(converted_cnf),
        ]
        w2u_proc = run_cmd(w2u_cmd, args.timeout)
        if w2u_proc.returncode != 0:
            detail = (w2u_proc.stdout + "\n" + w2u_proc.stderr).strip()
            return RunResult(instance, ttc_wmc, None, None, "W2U_FAIL", detail)
        divisor = parse_w2u_divisor(w2u_proc.stdout)
        normalize_projection_comments_for_approxmc(converted_cnf)
        multiplier = parse_w2u_multiplier(converted_cnf)

        approx_cmd = [
            str(args.approxmc),
            "--input",
            str(converted_cnf),
            "-v",
            "0",
            "-s",
            str(args.seed),
        ]
        approx_proc = run_cmd(approx_cmd, args.timeout)
        if approx_proc.returncode != 0:
            detail = (approx_proc.stdout + "\n" + approx_proc.stderr).strip()
            return RunResult(instance, ttc_wmc, None, None, "APPROXMC_FAIL", detail)
        unweighted_count = parse_approxmc_count(approx_proc.stdout)
        approx_wmc = unweighted_count * multiplier / (Decimal(2) ** divisor)
        rel = relative_error(ttc_wmc, approx_wmc)
        status = "OK" if rel <= args.tolerance else "MISMATCH"
        detail = f"smt={smt_path} w2u={converted_cnf}"
        return RunResult(instance, ttc_wmc, approx_wmc, rel, status, detail)
    except subprocess.TimeoutExpired as ex:
        return RunResult(instance, None, None, None, "TIMEOUT", str(ex))
    except Exception as ex:
        return RunResult(instance, None, None, None, "ERROR", str(ex))


def discover_instances(root: Path) -> list[Path]:
    suffixes = {".cnf", ".wcnf", ".gz", ".xz"}
    paths: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix in suffixes or path.name.endswith((".cnf.gz", ".wcnf.gz", ".cnf.xz", ".wcnf.xz")):
            paths.append(path)
    return sorted(paths)


def print_result(result: RunResult) -> None:
    ttc = "" if result.ttc_wmc is None else str(result.ttc_wmc)
    approx = "" if result.approxmc_wmc is None else str(result.approxmc_wmc)
    rel = "" if result.rel_error is None else str(result.rel_error)
    print(f"{result.status}\t{result.instance}\t{ttc}\t{approx}\t{rel}\t{result.detail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", help="Local Track4_PWMC directory. If omitted, clone the repo into --work-dir.")
    parser.add_argument("--repo-url", default=DEFAULT_REPO_URL)
    parser.add_argument("--track-dir", default=DEFAULT_TRACK_DIR)
    parser.add_argument("--work-dir", type=Path, default=Path("/tmp/ttc_track4_pwmc"))
    parser.add_argument("--ttc", type=Path, default=Path("build-static/ttc"))
    parser.add_argument("--w2u", type=Path, default=DEFAULT_W2U)
    parser.add_argument("--approxmc", default=None)
    parser.add_argument("--precision", type=int, default=7)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--max-files", type=int, default=10)
    parser.add_argument("--tolerance", type=Decimal, default=Decimal("0.50"))
    parser.add_argument("--pattern", default=None, help="Substring filter on instance path.")
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.ttc = args.ttc.resolve()
    if not args.ttc.exists():
        raise SystemExit(f"ttc binary not found: {args.ttc}")
    if not args.w2u.exists():
        raise SystemExit(f"weighted_to_unweighted.py not found: {args.w2u}")

    approxmc = resolve_executable(args.approxmc, APPROXMC_FALLBACKS)
    if approxmc is None:
        found = shutil.which("approxmc")
        approxmc = Path(found) if found else None
    if approxmc is None:
        raise SystemExit("ApproxMC not found. Pass --approxmc /path/to/approxmc.")
    args.approxmc = approxmc

    instances_root = ensure_instances(args, args.work_dir)
    instances = discover_instances(instances_root)
    if args.pattern:
        instances = [p for p in instances if args.pattern in str(p)]
    if args.max_files > 0:
        instances = instances[: args.max_files]
    if not instances:
        raise SystemExit(f"no instances found under {instances_root}")

    print("status\tinstance\tttc_wmc\tapproxmc_wmc\trel_error\tdetail")
    start = time.time()
    failures = 0
    for instance in instances:
        result = run_one(instance, args, args.work_dir)
        print_result(result)
        if result.status not in {"OK"}:
            failures += 1

    elapsed = time.time() - start
    print(f"# checked={len(instances)} failures={failures} elapsed={elapsed:.2f}s")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
