# Agent TODO: volesti single-polytope volume (and then sampling) in ttc

## Plan
- [x] 1. Understand the volesti oracle ttc uses for one polytope
- [x] 2. Build a standalone test harness (`tests/volesti/`) and reproduce the long-box bug
- [ ] 3. Fan out 10 test agents, one per area (done: e2e, annealing; 8 still running)
- [x] 4. Read the theory ([CV-PVA], [CV-GC], MATLAB) → `ttc-volume-vs-theory.md`
- [x] 5. Keep all testing on ≤ 60% of the CPU (cores 0–11, nice 19; a watchdog enforces it)
- [ ] 6. Consolidate a root-cause diagnosis from the agent reports
- [ ] 7. Fix: scale normalisation, rounding, a safe annealing step, scale-free volesti constants
- [ ] 8. Validate: unit suite + e2e suite, rebuild ttc, confirm no regressions on well-rounded bodies
- [ ] 9. **Sampling phase**: on the same polytopes, test that ttc's samples are uniform over each body
- [ ] 10. Commit after every step

---

## Step 1: Understand the oracle
**Did:** read `src/volume/lra_volume.cpp` and the volesti headers it includes.
**Found:**
- ttc calls `volume_cooling_gaussians<GaussianBallWalk, mt19937 seed 3>(P, ε/12, 10+n²)`
  (`lra_volume.cpp:1255`). This is Cousins–Vempala Gaussian cooling with a Gaussian ball walk.
- The polytope is not rounded first. The inner ball is pre-set by `reliableInnerBall`.
- The seed is fixed, so `--seed` does not change the volume estimate.
- `update_delta` takes the walk by value, so it has no effect. That doesn't matter, because
  the walk's constructor already computes the same δ.

## Step 2: Harness + reproducer (commit `1c35e90`)
**Did:**
- Wrote `tests/volesti/common.hpp`: body builders with exact volumes, plus the oracle run in
  exactly ttc's configuration.
- Wrote a doctest `Makefile`, and `test_smoke.cpp`.

**Found:**
- Cubes in 2–6D are within 3%.
- **Box [0,2]×[0,1000]: the oracle gives 7.78 (exact 2000).** ttc prints `s volume 7.85`.

## Step 3: Test agents
### e2e (commit `61eeb0d`)
**Did:** wrote 68 single-polytope `.smt2` files with exact volumes, and
`tests/volesti/run_smt2.py`.
**Found:** 35 of 67 cases fail.
- **Long boxes:** the estimate is far too low, with an error that grows with aspect ratio.
  Aspect 10 already fails in 5D (0.36×).
- **Side 1e-3:** the estimate is too low. ttc's clip lifts it to the lower bound.
- **Side 1e-6:** the run **hangs**. Cause: the absolute tolerance `tol=1e-7` in the
  `get_first_gaussian` bisection.
- **Side 1e6:** in 3D and 5D the oracle returns **inf** and ttc exits with code 1.
- **Equality constraints:** the volume is correctly 0.

### annealing (`tests/volesti/test_annealing.cpp`)
**Did:** wrote a bit-exact traced copy of the oracle, and computed exact per-phase ratios
(erf products for boxes).
**Found:**
- On the long box, all of the loss is in the **last phase** (a: 0.517 → 0). The estimate there
  is 2.33 against an exact 563, a factor of 0.0041. Phases 0–2 are only about 2% low each.
- **Cause A:** `get_next_gaussian` (`volume_cooling_gaussians.hpp:153-176`) tests the variance
  only on samples near the centre. So k keeps doubling, and the step jumps about 32 halvings
  at once, a jump whose exact relvar is infinite.
- **Cause B:** the ratio 1−1/n = 0.5 in 2D is never safe. The paper's Lemma 3 needs n ≥ 4.
- The "replace aᵢ by 0" branch is *not* the culprit.

**Validated patch:** cap each step at a′ ≥ a(1−α_max(n,C)), derived from the log-concavity of
aⁿZ(a). With it:
- long box 7.78 → 842
- cubes still within 3.5%
- cost: about 5–8× more phases.

The remaining long-box gap comes from ball-walk mixing:

| Walk (with cap) | Walk length | Long-box volume |
|---|---|---|
| Ball | 1400 | 1754 |
| CDHR | 1 | 1930 |

So **rounding or a better walk is still required**.

### gaussianwalk (`tests/volesti/test_gaussian_walk.cpp`, 8 cases, ~4 s)
**Did:** tested GaussianBallWalk as an MCMC sampler:
- uniformity of the sphere proposal
- the `is_in` convention
- δ
- stationary marginals against analytic truncated normals
- reversibility
- autocorrelation
- exact per-step relvar on the oracle's own schedule

**Found:** **the walk is correct.** All the stationarity and reversibility tests pass. It just
cannot mix on the long box:
- The Chebyshev centre is picked at the strip end (1, 999).
- δ = 2.83, so one end-to-end traversal takes about 9.5e5 steps. A schedule phase only reaches
  y ≈ −216.
- Lag-1 autocorrelation is 0.9989.
- Exact relvar of the oracle's schedule on the long box: 10^223944.
- A 10× longer walk still gives −86%. CDHR alone still gives −26 to −61%.
- `update_delta` taking the walk by value is cosmetic: results are bit-identical with it fixed.

**Conclusion:** the fix must be rounding and/or schedule repair, not a different walk.

## Step 4: Theory
**Did:** read [CV-PVA] (§2–4.2), [CV-GC] (Algorithm figure, assumptions), and the MATLAB
`Volume.m`, `ConvexBody.m` and `preprocess.m`. Wrote `ttc-volume-vs-theory.md`.
**Found:**
- volesti is a near line-by-line port of the MATLAB code.
- **ttc skips the mandatory rounding step.** It is [CV-PVA] Fig. 2 step 1. The MATLAB code's own
  usage example rounds `long_box` first.
- volesti adds absolute constants that the MATLAB code doesn't have: a bisection tolerance, and
  the `max(1,a)` term in the ball-walk step. These break scale invariance.
- Even the reference implementation is biased low and exceeds ε on simplices, so the oracle
  is heuristic.

## Step 5: CPU cap
**Did:** pinned the Claude process tree to cores 0–11 (12 of 20) at nice 19. Some processes
were started detached or reset their affinity, so a watchdog
(`scratchpad/cpu_guard.sh`) re-pins every 2 s any process under the Claude tree or with its
cwd/exe in the repo or scratchpad. Told all agents to run one process at a time.
