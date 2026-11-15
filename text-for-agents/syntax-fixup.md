# Syntax fix-up plan

*2026-11-15. Fixes for the SMT-LIB syntax failures found by the fuzzer (`fuzz-report.md`, §1).
Scope: input handling of QF_LRA volume counting. The algorithm and the volume engine do not change.*

## Plan
- [ ] S1. Catch parse errors and print a clean error instead of aborting (class B)
- [ ] S2. Inline `define-fun` definitions; count only free constants (class A, the big one)
- [ ] S3. Reject `define-fun-rec` / `define-funs-rec` cleanly
- [ ] S4. Decision: `to_real` under QF_LRA (strict error, or relax to QF_LIRA)
- [ ] S5. Decision: missing `set-logic` (cvc5's ALL, or default to QF_LRA)
- [ ] S6. Take the logic from cvc5's parsed state, not from a raw text search
- [ ] S7. Tests: `tests/volesti/syntax_smt2/` with exact volumes; rerun the fuzzer with the same
  seeds
- [ ] S8. Regression gate and commits (one commit per step)

Class C (timeouts on formulas with many overlapping polytopes) is performance work and stays out
of scope here.

---

## S1. Parse errors abort (class B, 89 CRASH)
**Problem:**
- `parser.parseFormula(...)` is called outside any try block at `src/main.cpp:1868`, and also at
  :405 and :477, the other engine entry points.
- A `cvc5::parser::ParserException` or `CVC5ApiException` therefore ends in `terminate called …`,
  exit code 134. Examples: `to_real` under `(set-logic QF_LRA)`, a type mismatch, an undeclared
  symbol.

**Fix:**
- Wrap every `parseFormula` call.
- On an exception, print `Error: parse error: <what()>` and return 1, the same exit code as the
  other `Error:` paths. (Check: `volume computation failed` returns 1.)

**Risk:** none. Behaviour changes only on inputs that used to abort.

## S2. `define-fun` (class A: 150 WRONG, about 556 UNSUPPORTED, 39 CRASH)
**Problem:**
- cvc5 records every `define-fun` as an assertion `(= f body)`, or `(= f (lambda (params) body))`,
  where f is a fresh constant (`SolverEngine::defineFunction` → `addDefineFunDefinition`).
- `getAssertions()` returns these assertions. ttc conjoins them into the formula
  (`src/parser.cpp:505`), and `collect` (`:623`) registers f and the lambda's bound parameters as
  counting variables.
- This happens **even when f is never used**:
  - 0-ary Real → an extra dimension pinned by `=` → `s volume 0`;
  - 0-ary Bool → a Boolean variable → the PACT engine (THEORY_BV error, or `s mc 1`);
  - with arguments → "unsupported literal", or THEORY_UF.

**Fix** (draft in `tests/volesti/fuzz_lra/suggested_fixes.diff`):
1. **Identify definitions robustly.** In the command loop (`parser.cpp` ~495), record the symbol
   of every command whose name is `define-fun`, from `cmd.getCommandName()` or by parsing
   `cmd.toString()`.
   - Do not rely on "not in `getDeclaredTerms()`", because the symbol manager may list defined
     symbols too. **Verify** which behaviour cvc5 has, with a 3-line input.
   - Mind `push`/`pop` scoping: a definition inside a popped scope must not apply.
     `getAssertions()` after the pop already reflects the scope, so match assertion
     `(= f body)` with f's symbol in the recorded set.
2. Split `getAssertions()` into definitions `(= f body)` and real assertions.
3. Substitute definitions into the real assertions until nothing changes. A body may use earlier
   definitions; definitions are non-recursive, so this terminates. Then
   `d_solver.simplify(...)` beta-reduces `(lambda …)` applications (the UF rewriter). Check this
   on `(define-fun g ((u Real)) Real (* 2 u))`.
4. Set `d_assertions` and `rawFormula` to the substituted list. That list must also feed `collect`
   (which runs on the unsimplified `rawFormula`) and the projection / weight lookups.
5. Change `collect`'s leaf test to also require `t.getKind() == cvc5::Kind::CONSTANT`, so bound
   `VARIABLE`s (lambda parameters) are never counted.

**Risk:**
- Every input containing a `define-fun` changes behaviour, which is the intent.
- Inputs without one must be **bit-identical** (gate in S8).

## S3. Recursive definitions
`define-fun-rec` / `define-funs-rec` cannot be inlined. Detect them in the same command loop and
fail with `Error: recursive function definitions (define-fun-rec) are not supported`, instead of
silently miscounting.

## S4. `to_real` under QF_LRA (63 aborts): decision needed
`to_real` belongs to the Reals_Ints theory, so cvc5 is right to reject it under QF_LRA. After S1
this becomes a clean error. Options:
- **(a) Strict** (default once S1 is done): the clean error plus a hint, "use (set-logic QF_LIRA)
  or drop to_real".
- **(b) Lenient:** on that specific `ParserException`, retry with the input's QF_LRA replaced by
  QF_LIRA, and print a warning.
  - The engine dispatch must still pick the LRA volume engine. Today it reads the text token, so
    S6 matters here.
  - Integer *variables* must still be rejected, because volume needs Real-only. Check
    `d_intVars` is empty.

**Recommendation:** (b). It is cheap and matches what people write.

## S5. Missing `set-logic`: decision needed
- Without `set-logic`, cvc5 uses logic ALL: numerals are Int, so mixed Int/Real terms such as a
  `define-fun` with Real parameters applied to `0` can fail to type-check.
- `extractLogicUpper` returns "", and dispatch falls back to its default.
- In the fuzz run, 234/258 such formulas were already OK. The failures were `define-fun` (fixed by
  S2) and 4 type mismatches.
- **Options:**
  - (a) Keep ALL. S1 turns the type mismatches into clean errors.
  - (b) If there is no `set-logic` and every declared symbol is Real, prepend `(set-logic QF_LRA)`
    before parsing.
- **Recommendation:** (b), behind a log line, "no set-logic: assuming QF_LRA".

## S6. Logic detection
`extractLogicUpper` (`main.cpp:1059`), and the similar search at :920 / :973 / :1002, take the
first textual `(set-logic` in the file. A `(set-logic QF_BV)` inside a comment or a string literal
would be picked up. The fuzzer's tricky comments happened not to trigger it.

**Fix:**
- Strip comments and string literals before searching, or better, read the logic cvc5 parsed
  (`d_solver.getLogic()` once `isLogicSet()`).
- The dispatch runs before parsing in some paths. If so, a tokenizer-based scan (skip `;`
  comments, `"..."` strings and `|...|` symbols) is the minimal fix.

## S7. Tests
**New `tests/volesti/syntax_smt2/`**: one file per construct, each with an `; exact-volume:`
header, run by `run_smt2.py` (`--dir` or `--filter`).

| File | Expected |
|---|---|
| `defun_real0_unused.smt2` | volume of the box (the constant is ignored) |
| `defun_real0_used.smt2` | `(define-fun c () Real 2.5)` used in a bound |
| `defun_bool0_unused.smt2` | volume, not `s mc` |
| `defun_bool0_used.smt2` | the Bool macro used as a predicate |
| `defun_affine.smt2` | `(define-fun g ((u Real) (v Real)) Real (+ u (* 2 v)))` |
| `defun_pred.smt2` | `(define-fun inb ((u Real) (lo Real) (hi Real)) Bool (and (<= lo u) (<= u hi)))` |
| `defun_nested.smt2` | a definition using another definition |
| `defun_push_pop.smt2` | a definition inside a popped scope must not apply |
| `defun_no_logic.smt2` | `define-fun` without `set-logic` |
| `defun_rec.smt2` | report-only: clean error (S3) |
| `to_real_qflra.smt2` | volume with S4(b), or a clean error with S4(a) |
| `int_literal_real_param.smt2` | the S5 behaviour |
| `logic_in_comment.smt2` | `; (set-logic QF_BV)` in a comment, real logic QF_LRA → volume (S6) |

**Fuzzer rerun**, same corpus (`--seed 1000..1003`, 4 × 750), with acceptance:
- WRONG = 0 and CRASH = 0.
- UNSUPPORTED only from `define-fun-rec`, plus the choices made in S4/S5.
- OK ≥ 99% of the non-timeout formulas.

Add the minimised reproducers in `tests/volesti/fuzz_lra/failures/` to the suite as fixed
regressions, with exact volumes from their headers' Monte-Carlo values or by hand.

## S8. Regression gate
After each step, all of the following must hold:
- `run_smt2.py`: 67/67, plus the new `syntax_smt2/` files.
- `run_union_e2e.py --seeds 10`: 38/38.
- `run_union_decomp.py`: 21/21.
- The default unit suite passes.
- **Bit-identity** for inputs without `define-fun`: the 16-run baseline (4 inputs × default/
  `--nogmp`/`--fullgmp`/`--sampler cdhr`) gives an unchanged `s volume`.

Commit one change per step, with the tests in the same commit. Record each step in this file, as
a **Did / Found** section under the plan.

## Order and effort
S1 (minutes) → S2 + S3 (about 1–2 h with tests) → S6 → S4/S5 once decided → S7 fuzz rerun →
S8 gate. Only S4 and S5 need your input.
