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
  - [ ] u-robust: T2 error injection (volume ±ε′, sampler TV η) → `test_union_robust.cpp`
  - [x] u-params: T3 lattice count vs volume + T4 Thresh/ε′/δ′/clipping →
    `test_union_lattice.cpp`, `test_union_params.cpp`
  - [ ] u-decomp: T6 decomposition soundness/completeness via `--dump-ine` →
    `tests/volesti/union_decomp/`, `run_union_decomp.py`
  - [ ] u-e2e: T7 SMT formulas with exactly known union volumes, 20 seeds, all modes, current vs
    start binary → `tests/volesti/union_smt2/`, `run_union_e2e.py`
- [ ] U5. Consolidate the diagnosis; fix what is wrong (with the user's go-ahead on cost-affecting
  changes)
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

