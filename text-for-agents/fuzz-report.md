# ttc fuzzing report

*2026-11-14/15. Overnight campaign `night` (23:42–04:37, binary `3fe72b1`) and morning campaign
`night2` (06:48–08:50, binary `31d7eff`). Each campaign ran at most 10 parallel jobs, all at
nice 19. Details and per-step logs are in `agents-todo-union-e2e.md`.*

## Bottom line
- **Volume accuracy holds up.** 3600 runs on random Boolean unions with exactly known volume, in
  every GMP mode, gave **0 accuracy failures**. The worst error was 6.5% against ε = 0.8.
- **The fuzzing found four bugs:**
  1. **`define-fun` breaks ttc**, even when the defined symbol is never used: wrong `s volume 0`,
     a switch to the wrong engine, or errors. **Open.**
  2. **Parse errors abort the process** (signal 6) instead of printing an error. **Open.**
  3. **The default volume oracle could hang** (simplices, `--seed 10`). **Fixed** in `31d7eff`.
  4. A harness argument bug kept 12 union batches from running in night 1. Fixed; they ran in
     night2.
- **Without** `define-fun`, `to_real` under QF_LRA, or a missing `set-logic`, the syntax fuzz
  gives **1992/2007 OK and 15 timeouts, with 0 wrong answers and 0 crashes.**

## 1. Syntax fuzz: 3000 tiny satisfiable QF_LRA formulas
**Setup:**
- Tools: `tests/volesti/run_fuzz_lra.py` with generator `tests/volesti/fuzz_lra/gen_fuzz.py`.
- Formulas: 1–4 Real variables, each inside a bounding box, with randomly chosen syntax features.
- Each formula is kept only if a Monte-Carlo evaluator confirms it is satisfiable with positive
  volume.
- ttc's answer is checked against that Monte-Carlo volume.
- Run: night2, 4 × 750 formulas (seeds 1000–1003), 20 s timeout per formula.

| Class | Count | Meaning |
|---|---|---|
| OK | 2187 (73%) | volume agrees with Monte Carlo |
| UNSUPPORTED | 556 | ttc printed an explicit error |
| WRONG | 150 | a finite answer far from Monte Carlo (all of them `s volume 0`) |
| CRASH | 89 | abort, a missing `s volume`, or the wrong engine |
| TIMEOUT | 18 | over 20 s |

### Failure classes
| # | Trigger | Effect | Count | Root cause | Fix (status) |
|---|---|---|---|---|---|
| A | `define-fun`, **even if unused**: 0-ary Real | `s volume 0` | WRONG 150; 150/200 of the formulas with this feature | cvc5 stores every `define-fun` as an **assertion** `(= f body)`, with f a fresh constant (`SolverEngine::defineFunction`). ttc conjoins `getAssertions()` into the formula (`src/parser.cpp` ~505), and `collect` registers f as a counting variable. f becomes an extra dimension pinned by `=`, so the volume is 0 | Remove the definitions from the assertions and substitute them; count only free `CONSTANT`s (**open**, diff in `fuzz_lra/suggested_fixes.diff`) |
| A | 0-ary Bool `define-fun` | routed to the PACT engine: `logic … doesn't include THEORY_BV`, or `s mc 1` without `set-logic` | 231 + 28 (CRASH) | the same cause; f is counted as a Boolean variable | same fix (**open**) |
| A | `define-fun` with arguments | `unsupported literal`, or THEORY_UF | ~370 | the same cause; the lambda's bound parameters carry symbols and are also collected | same fix (**open**) |
| B | `to_real` under `(set-logic QF_LRA)` | `terminate called … ParserException`, rc −6 | 63/65 | `to_real` belongs to Reals_Ints, not QF_LRA, so cvc5 is right to reject it. But `parser.parseFormula` (`src/main.cpp` ~1868) is outside any try block | Catch it and print `Error: parse error: …` (**open**). **Decision for you:** whether to accept `to_real` by switching to QF_LIRA |
| B | no `set-logic` + Int literal passed to a Real parameter | rc −6, "argument type is not …" | 4 | logic ALL types numerals as Int | same catch (**open**). **Decision for you:** whether to treat a missing `set-logic` as QF_LRA |
| B | undeclared symbol (a shrinking artifact) | rc −6 | 2 | same abort path | same catch (**open**) |
| C | heavy Boolean/`ite` structure | timeout | 18 | 12–78 overlapping polytopes; a performance issue, not correctness | — |

### What works
Leaving out formulas that use `define-fun`, `to_real` under QF_LRA, or no `set-logic`, **1992/2007
are OK and 15 time out.** Every remaining feature (≥ 20 uses each) gives no wrong answer, crash or
unsupported error:

| Area | Features |
|---|---|
| Declarations | `declare-fun`, `declare-const`, quoted and special symbols (`\|tab q\|`, `\|p"q\|`) |
| Connectives | `and`, `or`, `not` (nested), `=>`, `xor`, Bool `ite`, Bool `=` |
| Real `ite` | inside arithmetic (lifted since `23befd1`) |
| `let` | parallel and shadowing, Bool and Real |
| Comparisons | comparison chains (`< a b c`), n-ary and negated `=`, `distinct` (2 and n arguments) |
| Arithmetic | unary and n-ary `-`, n-ary `+`, `*` by a constant on either side, `/` by a constant, rationals, decimals, large and small constants |
| Atoms | constant atoms, duplicate atoms and asserts |
| Commands | `:named` (including references), `push`/`pop` (including unclosed), `set-info`/`set-option`, `check-sat` in the middle of the file, `get-model`/`get-value`, `exit` |
| Lexical | tricky comments, string literals in `set-info` |

## 2. Exact-volume union fuzz: 3600 ttc runs
**Setup:**
- Tools: `tests/volesti/fuzz_union/gen_bool_unions.py` and `run_batch.py`.
- Formulas: random `and`/`or`/`not`/`xor`/`=>`/`ite` combinations of axis-aligned atoms inside a
  box. The exact volume is the sum over compressed-grid cells, computed with fractions.
- Half the formulas are written in sheared coordinates x = M·y (integer M, volume × |det M|), so
  their atoms are not axis-aligned.
- Dims 2–6; sides 1e-3 to 1e3.
- 25 formulas per batch, 3 ttc seeds each, 240 s timeout.

| Binary | Variant | Dims | Runs | OK | Accuracy failures | Mean err | sd | Worst \|err\| | Median time |
|---|---|---|---|---|---|---|---|---|---|
| `3fe72b1` | default | 2–6 | 1800 | 1769 | 0 | +0.3 … +0.8% | 1.3–1.6% | 6.5% (5D) | 1–19 s |
| `3fe72b1` | `--sampler cdhr` | 3–4 | 450 | 450 | 0 | +0.2% | 1.5% | 4.8% | 4 s |
| `3fe72b1` | `--epsilon 0.3 --delta 0.1` | 2–3 | 450 | 445 | 0 (vs 0.3) | +0.1 … +0.3% | 0.7% | 2.3% | 8–27 s |
| `31d7eff` | `--nogmp` | 2–4 | 450 | 450 | 0 | +0.7% | 1.4% | 5.3% | 1.6 s |
| `31d7eff` | `--fullgmp` | 2–3 | 450 | 438 | 0 | +0.4% | 1.3% | 4.1% | 23.5 s |

- **Not OK = timeouts only** (240 s): 31 in the 4–6D default runs, 5 at ε=0.3 in 3D, and 12 with
  `--fullgmp`.
- `--fullgmp` is about 15× slower than the default.
- The small positive mean is the volume oracle's own bias, within its ε/12 budget. It is the same
  on single boxes; see the union e2e results.

## 3. Seeded end-to-end sweeps
- **Single-polytope suite** (`run_smt2.py --seeds 20`, 67 formulas with exact volume):
  - night 1 (`3fe72b1`): 1334/1340; **6 TIMEOUTs, all `--seed 10` on 2D/3D simplices**.
  - After the fix (`31d7eff`): **1340/1340**.
  - Cause: three bugs in volesti's CDHR oracle walk, fixed in `31d7eff`:
    - a facet through the current point was ignored, so the chain left the polytope;
    - the incremental slack updates drifted;
    - the Gaussian chord samplers had unbounded rejection loops on far-tail chords.
- **Union suite** (`run_union_e2e.py`, 38 formulas with exact volume):
  - 10 seeds: 38/38, 0 failing runs, |mean err| 0.3–1.0%.
  - 100 seeds (night2, stopped early): 29 formulas finished at ε=0.8 and 4 at ε=0.3, all PASS.

## 4. Timing, before vs now (side by side, 10 random syntheticLRA instances)
Before = `update` @ `75e05e0`; now = `31d7eff`. Both binaries ran concurrently on each instance.

| | Before | Now |
|---|---|---|
| solved within 900 s | 4/10 | 7/10 |
| time, when both solved | — | 0.50–0.64× before (the trivial dim-3 instance: 1.08×) |
| total wall time (timeout = 900 s) | 6258 s | 4864 s |
| volume, when both solved | — | 0.973–1.011× before |

**No slowdown: the current build is 1.5–2× faster.** Night 1's 5–15× gap against the archived
`ttc_correct` times comes from the archive being cluster runs, not from ttc.

## Reproducers
Minimised inputs are in `tests/volesti/fuzz_lra/failures/`. Each one's header records the class,
the signature, the Monte-Carlo volume and ttc's output:

| Failure | File |
|---|---|
| A, Real `define-fun` | `WRONG_volume_zero_vs_MC_ratio_N__0409.smt2`: an unused `(define-fun K1 () Real (- 4240000))` → `s volume 0` instead of about 1.1e6 |
| A, Bool `define-fun` | `UNSUPPORTED_Error_approximate_counting_failed_*.smt2`, and `CRASH_rc_N_no_s_volume_s_mc_N__0082.smt2` (no `set-logic` → `s mc 1`) |
| A, with arguments | `UNSUPPORTED_Error_volume_computation_failed_unsupported_lite__0268.smt2` |
| B, parse-error abort | `CRASH_rc_N_Symbol_to_real_not_declared_as_a_variable__*.smt2`, `CRASH_rc_N_argument_type_is_not_*.smt2` |
| C, timeout | `TIMEOUT_timeout_polytopes_N__*.smt2` |

## How to rerun
```sh
python3 tests/volesti/run_fuzz_lra.py --ttc build-gcc/ttc --count 750 --seed 1000 --jobs 2
python3 tests/volesti/fuzz_union/run_batch.py --ttc build-gcc/ttc --work /tmp/bu --seed 100 \
        --count 25 --seeds 3 --dims 2,3,4 --extra-args='--nogmp'
python3 tests/volesti/night_runner.py --root DIR      # overnight queue (DIR/jobs/*.json)
python3 tests/volesti/night_summary.py --root DIR     # report
```

## Open items
1. **Fix A (`define-fun`)** and **fix B (catch parse errors)**. Both are small; the drafts are in
   `tests/volesti/fuzz_lra/suggested_fixes.diff`. Not applied yet, because applying them needs a
   build and a test rerun.
2. **Decisions for you:**
   - accept `to_real` under QF_LRA by switching to QF_LIRA;
   - treat a missing `set-logic` as QF_LRA.
3. **Performance:** formulas with many overlapping polytopes (class C), and `--fullgmp` at about
   15× slower.
