# Agent TODO: union accuracy (KR25 Alg. 1, end to end)

Test plan: `text-for-agents/ttc-union-test-plan.md` (components C1–C8, tests T1–T7).
Manuscript: `text-for-agents/kr25-ssm/sections/{algo,analysis,proof-appendix}.tex`.
Single-polytope volume and sampling work: `agent-todo.md`.

## Plan
- [x] U1. Read the manuscript; write the component test plan (`64fb1d8`)
- [x] U2. Bit-identical refactor: union loop + store + Thresh/GetPrecision →
  `src/volume/union_estimator.{hpp,cpp}`, result structs → `volume_result.hpp` (`200276c`)
- [x] U3. Test harness links `union_estimator.o` + GMP (`tests/volesti/Makefile`)
- [~] U4. Fan out 5 agents (CPU cap: 12 cores, nice 19, ≤ 2 processes each) — all 5 launched:
  - [ ] u-core: T1 estimator with exact components + T5 store mechanics → `test_union_core.cpp`
  - [x] u-robust: T2 error injection (volume ±ε′, sampler TV η) → `test_union_robust.cpp`
  - [x] u-params: T3 lattice count vs volume + T4 Thresh/ε′/δ′/clipping →
    `test_union_lattice.cpp`, `test_union_params.cpp`
  - [x] u-decomp: T6 decomposition soundness/completeness via `--dump-ine` →
    `tests/volesti/union_decomp/`, `run_union_decomp.py`
  - [ ] u-e2e: T7 SMT formulas with exactly known union volumes, 20 seeds, all modes, current vs
    start binary → `tests/volesti/union_smt2/`, `run_union_e2e.py`
- [~] U5. Consolidate the diagnosis; fix what is wrong. Changes the user approved on 2026-11-14:
  - [x] U5b. GetPrecision 4r/η → 8r/η in code (`01ac059`) and in the paper (kr25-ssm algo/analysis/
    proof-appendix: edited, **not committed**, because the paper repo has the user's own uncommitted
    edits)
  - [x] U5c. Volume-oracle seed derived from `--seed` and the polytope index (`92395cf`)
  - [x] U5d. Variable names in the `--dump-ine` header (`3fe72b1`)
  - [ ] U5e. Lift a Real-valued `ite` inside arithmetic atoms
- [~] U8. New requests: tiny-formula syntax fuzzing (agent u-fuzz), general speed-ups (agent
  u-perf, own git worktree)
- [~] U9. Overnight campaign (≤ 10 jobs), started 23:42; read it in the morning with
  `python3 tests/volesti/night_summary.py --root <night dir>`
- [ ] U6. Validate (all suites, e2e single + union, benchmark), report in
  `ttc-union-vs-theory.md`
- [ ] U7. Commit after every step

---

## U1. Read the manuscript, write the test plan (`64fb1d8`)
**Did:** read `algo.tex` (Alg. 1 `ttc`, Alg. 2 `GetPrecision`) and `analysis.tex` (lemmas
`decomposition`, `process`, `integerlattice`, `nonconv`, `lattice`, `numlattice`). Checked the
code against both algorithms.

**Found:**
- `lraSampleThreshold`, ε′ = ε/12, δ′ = δ/(4m), the halving test `p ≥ Thresh/t` (computed in log
  space) and the output |X|/(p·Thresh·10^{bn}) all match Alg. 1.
- `precisionFromCubes` matches Alg. 2: γ = 64mn/η·√log(4r/η), with η = ε/96.
- **Concern for T5:** `SampleStore::removeInside` tests `A·x − b ≤ 0` in floating point, with
  tolerance 0. Lattice points on a facet shared between polytopes can be misclassified, which would
  double-count overlaps. This is the same class of bug as the sampler's lattice-membership bug,
  fixed in `08026a2`.

## U2. Refactor (`200276c`)
**Did:**
- Moved `SampleStore`, `PolytopeVolume` and `runUnionAlgorithm` (a header template on the
  coordinate type and the generator) into `union_estimator.hpp`.
- Moved `lraSampleThreshold` and `precisionFromCubes` into `union_estimator.cpp`.
- Moved the result rows into `volume_result.hpp`.
- `lra_volume.cpp` includes the new header.

**Found:** bit-identical. The fresh 16-run baseline (4 inputs × default/`--nogmp`/`--fullgmp`/
`--sampler cdhr`) is unchanged.

## U3. Test harness
**Did:** `tests/volesti/Makefile` now compiles `union_estimator.o`, adds `-I src/ddnnf` for
`logger.hpp`, and links `-lgmpxx -lgmp`. `test_smoke` still passes.

## U4. Agents
**Did:** launched u-core (T1+T5), u-robust (T2), u-params (T3+T4), u-decomp (T6) and u-e2e (T7).
Each owns its own test files, makes no edits under `src/`, makes no commits, and runs at most 2
niced processes.

### u-params (T3 + T4): `test_union_lattice.cpp` (7 cases, 4.6 s), `test_union_params.cpp` (6 cases, 6.7 s)
**Did:**
- Recomputed Alg. 2 by hand over 11760 grid points.
- Checked the Azuma tail with natural log against log₁₀.
- Ran the proof steps of `lem:nonconv` numerically.
- Enumerated |Q∩L| exactly on 13 unions: L-shape, frame, staircase, triangles tiling a square or a
  hexagon, 3D L, cube∪simplex, …
- Checked Thresh against Alg. 1 line 5, and checked ttc's wiring of ε′, δ′ and b.
- Checked clipping soundness on 160 bodies, and scale extremes (dims up to 40, sides 1e±9).

**Found:**
- **All correct in code:**
  - `precisionFromCubes` matches Alg. 2 exactly. The natural log is the right one; with log₁₀ the
    Azuma bound fails.
  - Thresh matches the formula and satisfies both inequalities of the `numlattice` proof.
  - The clip bounds are sound: γ ≤ true inradius, R ≥ max vertex distance.
  - `lem:lattice` holds at ttc's b with |err| ≤ 2.3e-5 against η ≈ 3e-3.
  - At extreme scales ttc throws rather than returning a wrong value.
- **Proof-text gap in `lem:nonconv`** (analysis.tex ~l.193). The step (1+κ/γ)ⁿ ≤ (1+η/32mn)ⁿ cites
  √u+√v ≤ √(2(u+v)), which only gives κ ≤ 2√ln(8r/η), but Alg. 2 uses ln(4r/η). The code-level
  check fails at 2913 of 720k points.
  - The lemma's conclusion still holds because of slack elsewhere, so there is **no accuracy impact**.
  - **Fix:** use 8r/η in Alg. 2 and at `union_estimator.cpp:17`; validated in an overlay, with 0
    failures. **Decision for the user (paper + code).**
- **GetPrecision is conservative:** ttc's b is 2–4 digits above the smallest η-accurate b.
  - It costs almost nothing: C is tracked in logs, and the extra halvings run on an empty store.
  - Exception: with `--nogmp`, double coordinates resolve 10⁻ᵇ only while |x|·10ᵇ < 2⁵³.

### u-decomp (T6): `run_union_decomp.py` + `union_decomp/d01–d21.smt2` (`4529f24`)
**Did:**
- Dumped the polytopes with `--dump-ine` and checked F(x) ⇔ ∃i: x ∈ P_i using its own QF_LRA
  evaluator.
- Points: uniform in a bounding box, inside each P_i, inside each DNF region, and just across every
  facet. Interior mismatches were confirmed with cvc5.
- Inputs: 21 formulas plus 5 syntheticLRA instances.

**Found:**
- **No unsound polytope anywhere** (no P_i ⊄ sol(F)).
- **Overlap ratio** Σvol(P_i)/vol(Q) ranges from 1.0 to 1.48. That overlap is what the union loop
  has to handle.
- **Three completeness bugs, all pre-existing** (the start binary behaves the same):
  1. **`distinct` / ≠:** a negated equality threw, and `parser.cpp:805-807` swallowed the
     exception, so ttc printed **`s volume 0` with exit 0** (true volume 1).
  2. **Tiny bodies dropped:** cddlib's absolute tolerances (about 1e-7) reported full-dimensional
     polytopes of width about 1e-7 as flat, and they were dropped (volume 0).
  3. **Real-valued `ite`** inside an arithmetic atom is unsupported.
- The `.ine` dump has no variable names, and its column order is the DFS first-occurrence order.

## U5a. Decomposition fixes (`33c8be4`)
**Did:**
- a ≠ b now adds no constraint (it removes only a hyperplane, which has measure zero).
- `canonicalizePolytope` runs cdd in the certified Chebyshev frame x = c + r·y with unit-norm rows,
  then maps back.
- `parser.cpp` rethrows if both decomposition attempts fail, instead of returning an empty cover.

**Found:**

| Input | Exact | Before | After |
|---|---|---|---|
| d09 `distinct` | 1 | 0 | 1.037 |
| d19 frame at 1e-7 scale | 1.2e-13 | 0 | 1.229e-13 |
| d20 thin union | 1.5e-7 | 0 | 1.535e-7 |
| d07 Real-valued `ite` | 2 | 0 (exit 0) | explicit error |

- `run_union_decomp`: **19/21 pass**. d07 fails by design (explicit error); d14 (free Bools) is
  report-only.
- `run_smt2`: 67/67, and the equality inputs still give 0.
- **Open, for the user:**
  - lifting a Real-valued `ite` into Boolean structure;
  - writing the variable names into the `.ine` header.

### u-robust (T2): `test_union_robust.cpp` (7 cases, 162 s; slow sweeps tagged)
**Did:**
- Test unions: 2D parallelotopes (pair90, chain8 with overlap ρ = 4.7, nested, rotated slabs3,
  crossX), plus slabsThin, chain5d and skinny5d (rotated, aspect 100) for ttc's components.
- Sampler: exactly uniform on the lattice points. All lattice counts enumerated exactly.
- Predicted mean E[est] = Σ V̂_i·Q_i(D_i), where D_i = P_i minus all later polytopes.

**Found:**
- **No union bug.** Every exact-component prediction holds (|z| ≤ 4).
  - Baseline bias ≤ 0.12%; sd 0.6–0.9% at ε=0.8.
- **Volume errors:** the bias is linear, Σe_i·|D_i|/|Q|, whatever the overlap pattern. The union
  tolerates e ≈ 5.9 ε′: e = 5ε′ gives 0/200 failures.
- **Sampler bias** (TV = total-variation distance from uniform):
  - |bias| ≤ ρ·TV. Runs pass at TV = ε′ and fail at 2ε′.
  - So the tolerated TV is about (ε/2)/ρ, six times the plan's ε′/m estimate.
  - The lemma's pointwise (1±ε′) assumption gives |b| ≤ ε′, as measured.
- **ttc's current components are inside the tolerance:**
  - oracle |e| ≤ 0.42 ε′ (ε=0.8);
  - sampler mass on each D_i within ±5%, union bias ≤ 0.6%;
  - full union: **0 failures** (100 seeds at ε=0.8, 20 at ε=0.3).
- **`--no-round` fails:** crossX −9.4%, skinny5d −28.6% (beyond ε/2). This is the old restart bias.
- **Design note:** the oracle's fixed seed 3 gives congruent polytopes identical errors (+0.0282).
  Volume errors are therefore fully correlated across polytopes and across ttc `--seed`s, so they
  never average out. Within budget here. **Decision for the user:** derive the oracle seed from
  `--seed`.

## U5b–d. Approved changes
- **U5b** (`01ac059`): `precisionFromCubes` now uses √ln(8r/η).
  - `test_union_lattice` is 7/7 and e2e 67/67.
  - Paper: algo.tex l.73/88, analysis.tex l.127/186/251, proof-appendix l.162/172 (cor:linf → 64r/ε).
    Edited, not committed.
- **U5c** (`92395cf`): per-polytope oracle streams, splitmix64(seed, index). Retries are
  attempt-salted. e2e at 3 seeds: 201/201.
- **U5d** (`3fe72b1`): `.ine` files carry `* variables: x y ...`. `run_union_decomp.py` reads the
  header.

## U9. Overnight campaign (started 2026-11-14 23:42, stops launching at 08:00)
**Where:**
- Runner: `tests/volesti/night_runner.py`, detached with setsid/nohup.
- Root: `<scratchpad>/night` (= `/tmp/claude-1000/-home-arijit-solvers-ttc-cpp/c5fd8853-63a5-4cce-a852-c234836b8395/scratchpad/night`).
- Frozen binary `night/bin/ttc` at `3fe72b1` (before Real-`ite` lifting). Frozen test binaries in
  `night/bin/tests/`.

**Workers:** 6 until 00:30, then 10 (`night/workers.txt` overrides). All jobs run at nice 19.

**126 jobs, round-robin interleaved:**

| Family | Jobs | Contents |
|---|---|---|
| `bu_*` | 48 | `fuzz_union/run_batch.py`: 25 random Boolean unions × 3 seeds each, with exact volume (grid cells, half sheared). Variants: dims 2–6; `--nogmp`, `--fullgmp`, `--sampler cdhr`, `--epsilon 0.3 --delta 0.1`. 3600 formula runs in total |
| `syn_*` | 60 | random syntheticLRA instances, current binary, `-v 2`, 1500 s limit; compared with the archived `ttc_correct` volume and time |
| `slow_*` | 17 | every `slow`-tagged unit sweep (`test_union_core` excluded: still being written) |
| `e2e_single_20seeds` | 1 | `run_smt2.py --seeds 20` |

Syntax-fuzz and union-e2e jobs are to be added when the u-fuzz and u-e2e agents deliver their
runners.

**Morning:** `python3 tests/volesti/night_summary.py --root <night dir>`.

