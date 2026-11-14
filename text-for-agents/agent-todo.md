# Agent TODO: volesti single-polytope volume (and then sampling) in ttc

## Plan
- [x] 1. Understand the volesti oracle ttc uses for one polytope
- [x] 2. Build a standalone test harness (`tests/volesti/`) and reproduce the long-box bug
- [x] 3. Fan out 10 test agents, one per area (all 10 reported)
- [x] 4. Read the theory ([CV-PVA], [CV-GC], MATLAB) → `ttc-volume-vs-theory.md`
- [x] 5. Keep all testing on ≤ 60% of the CPU (cores 0–11, nice 19; a watchdog enforces it)
- [x] 6. Consolidate a root-cause diagnosis from the agent reports
- [x] 7. Fix: scale normalisation (`6494cf2`), rounding + `--vol-oracle` (`882cf3f`), volesti rounding bugs
- [~] 8. Validate: e2e 67/67 pass; unit suite running; syntheticLRA regression check pending
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

### elongated (`tests/volesti/test_elongated.cpp`, 10 cases, ~12 s; slow sweeps tagged `slow`)
**Did:** aspect-ratio sweeps over boxes (one long axis, several long axes, a thin axis),
stretched simplices and rotated long boxes, in n = 2..10, plus walk-length scaling.

**Found — failure envelope:**
- One long axis fails from aspect 5 at n=2,3, and from about 10 at n=5..10.
- Stretched simplices fail from aspect 2–5.
- At aspect ≥ 100 it **always underestimates** (−0.67 to −1.0).
- In the 10–100 range, single runs can also overshoot (up to +9.1).

**Found — walk length** (R=100, 2D): mean error −0.63 at walk length 14, −0.20 at 140, −0.08 at
1400. That run took 513 s, and the required length grows like aspect².

**Found — side bug:** `max_inscribed_ellipsoid.hpp:132` uses an eigen-solver that throws
`std::invalid_argument` for n < 3.

**Validated patch:** map the maximum inscribed ellipsoid to the unit ball before cooling.
Iterate up to 5×, stopping once the axis ratio is ≤ 6, and multiply the result by |det|.
Guard against a non-finite transform or a centre outside P.
- All elongated cases: ≤ 0.06 mean error up to aspect 1e4 in n = 2..8.
- Small cubes (side 1e-3): −0.93 → −0.017.
- Runtime drops from 11.5 s to 3.6 s.

### altoracles (`tests/volesti/test_alt_oracles.cpp`, ~19 s; full sweep tagged `slow`, env ALT_*)
**Did:** compared volesti's volume oracles, 10 seeds each, with and without
max-inscribed-ellipsoid rounding (R+):
- CG (cooling Gaussians) with the Ball walk (ttc) and with CDHR
- CB (cooling balls) with the Billiard, Accelerated Billiard and CDHR walks
- SOB (sequence of balls)

**Found — mean relative error:**

| Body | CG-Ball (ttc) | R+CG-Ball | R+CB-Bill |
|---|---|---|---|
| [0,2]×[0,1000] | −0.87 | +0.002 | +0.007 |
| rot-box5-asp100 | −0.72 | +0.011 | +0.001 |

- **Without rounding no oracle is reliable** on rotated or skinny bodies. CDHR only looks fine on
  axis-aligned boxes; rotated, it gives −0.89.
- AcceleratedBilliard without rounding **overshoots by up to 1e13**.
- **Root cause:** every walk's step size scales with the inner radius r, not with the body size.
  - Ball walk: 4r/√(max(1,a)n).
  - Billiard walk: L = 2√n·r (`compute_diameter.hpp:39`).

**Found — speed:**
- After rounding, CB-Billiard is about 100–2000× faster than ttc's CG-Ball (cube20: 0.024 s vs
  47 s) with sd about 0.05–0.07.
- R+CG-CDHR: sd about 0.04, about 100× faster.

**Found — bug:** volesti's own rounding throws in **d=2** (Spectra `nev` must be ≤ n−1) at
`inscribed_ellipsoid_rounding.hpp:94` and `max_inscribed_ellipsoid.hpp:132`. Patched in an
overlay: use the Eigen eigensolver when d < 3 or when Spectra fails.

### invariance (`tests/volesti/test_invariance.cpp`, ~43 s; full size grid tagged `slow`)
**Did:**
- Scale s ∈ {1e-6..1e6}, translation |t| up to 1e9, rotations.
- Inner balls ×0.1..×10 and off-centre.
- Each oracle call runs in a forked child with a CPU cap, so hangs count as timeouts.

**Found:**
- **62 of 94 grid cells fail:**
  - s=1e-6: timeout
  - s=1e-3: −0.82 to −1
  - s=1e6: `inf`
- Translation and rotation are fine when the oracle gets an exact ball.
- An inner radius ×0.1 gives −0.947 on simplex5.

**Causes:**
- (1) The `get_first_gaussian` bisection is absolute: it starts at `upper=1` with `tol=1e-7`.
  - small bodies: it never terminates
  - large bodies: a₀ gets stuck near 1e-7
- (2) The step 4r/√(max(1,a)n) is not unit-consistent.
  - small s: the step shrinks like s²
  - large s: every proposal is rejected, `a_vals[0]=0`, and the result is `(π/0)^{n/2}=inf`
- volesti's `ComputeInnerBall` also has absolute tolerances. ttc's `reliableInnerBall` already
  guards against them.

**Validated patch:** rescale to unit inner radius inside `volume_cooling_gaussians`, then multiply
by rⁿ. Every scale then gives identical results, and only elongation failures remain (26/94).

### theory (notes: `scratchpad/theory/notes.md`)
**Found:**
- Reproduced 7.78395 bit-exactly.
- Volesti stops cooling at a=7.4e-3 instead of the a ≲ 1/(2R²) = 2e-6 the theory needs. The last
  ratio's true value is 48.6, the estimate is 3.01, and the true Var/E² is about 10^800.
- volesti's own `rounding_test.cpp` rounds skinny bodies first, and none of its volume routines
  rounds internally.
- `inscribed_ellipsoid_rounding` is a port of the MATLAB `preprocess.m`.
- Max-ellipsoid or SVD rounding + ttc's oracle gives errors within ±4% for n up to 30.

### shapes (`tests/volesti/test_standard_shapes.cpp`, ~70 s; slow set `-ts=slow --no-skip`)
**Did:** 20 seeds plus the deterministic ttc call on:
- cubes n=1..20
- simplices n=2..15
- cross-polytopes 2..8
- rotated cubes, products of simplices, regular polygons
- Birkhoff B₃
- boxes with redundant rows and with row scales of 1e±6

**Found:**
- Everything passes except **standard simplices from n=7**. The failure is **spread, not bias**:
  sd 0.104 at n=7 and 0.176 at n=15, against the 0.0667 target.
- No shape is consistently low.
- The cause is ball-walk correlation. A simplex has R/r ≈ n, and rounding cannot change that.
- **GaussianCDHRWalk fixes every simplex and is 2.5–8× faster.** With it, simplex15 goes from sd
  0.176 at 208 s to sd 0.064 at 27 s. In ttc this is `--vol-oracle cg-cdhr`.

### rounding (`tests/volesti/test_rounding.cpp`, 9 cases, ~26 s; slow method comparison tagged `slow`)
**Did:** round + ttc oracle on:
- the long box
- aspect 1..1e4 in d=2/5/10
- rotated long boxes and long simplices
- the slab [0,1e-4]×[0,1]^4
- cubes, to check rounding does not hurt them

**Found:**
- Max-ellipsoid rounding fixes every hard case (max error 0.05).
- On cubes rounding is the identity.
- The rounding itself costs ≤ 0.01 s and makes the oracle 3–5× faster.
- LOG_BARRIER, VOLUMETRIC and VAIDYA rounding are equally good.
- **SVD rounding is bad:** −0.28/−0.90 on simplices, and it **hangs** on the slab because the
  `svd_rounding.hpp` retry loop has no cap.
- Confirmed vol(P) = |det T|·vol(P_rounded).

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

## Step 7: Fix
**Did:**
1. `6494cf2` (volesti):
   - Rescale to unit inner radius inside `volume_cooling_gaussians` and multiply by rⁿ.
   - Relative a₀ bisection tolerance.
   - Spectra → dense eigensolver fallback for d<3.
2. `882cf3f` (ttc):
   - New `src/volume/polytope_volume_oracle.hpp`, containing `roundForVolume` (max-inscribed-ellipsoid
     rounding, falling back to the unrounded body on failure) and `runVolumeOracle`.
   - `lra_volume.cpp` phase 1b rounds a copy, re-certifies its ball with `reliableInnerBall`, and
     multiplies the oracle result by |det T|. The original body keeps its radii for clipping,
     GetPrecision and sampling.
   - New flags: `--vol-oracle {cg-ball (default), cg-cdhr, cb-billiard}` and `--no-round`.
   - volesti preprocessing fixes:
     - Catch Spectra's exception on round bodies.
     - Make the `Computations failed` errors actually throw.
     - Fix the swapped-arguments overload.
3. Tests: `common.hpp` gains `ttcVolume` / `ttcVolumeSeeded`, which run ttc's real pipeline.
   `ttcOracle` stays the raw volesti call.

**Found:**
- e2e: **67/67 pass** (was 32/67). Worst error 0.078. Wall time 22.5 s (was 76 s).
- Long box: 7.78 → 2046.8 (+2.3%).
- Side 1e-6: hang → +3%.
- Side 1e6: inf → +3%.
