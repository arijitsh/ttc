# How close is ttc's single-polytope volume to Cousins–Vempala?

*2026-11-14. Scope: the per-polytope volume oracle only (not the union algorithm).*

## Sources read
- **[CV-PVA]** Cousins & Vempala, *A practical volume algorithm*, Math. Prog. Comp. 2016
  (`Cousins-Vempala2016_Article_APracticalVolumeAlgorithm.pdf`). Fig. 2 is the volume
  algorithm, Fig. 3 rounding, §3.2 annealing schedule, §3.3 convergence, §4.2 accuracy.
- **[CV-GC]** Cousins & Vempala, *Gaussian cooling and O\*(n³) algorithms for volume and
  Gaussian volume* (`CV-gaussian-cooling.tar.gz` → `vol3_journal.tex`). Fig. "Volume(K,ε)" at
  tex line 206.
- **[MATLAB]** Cousins' reference code (`Volume-and-Sampling-master.zip`): `Volume.m`,
  `ConvexBody.m`, `preprocess.m`.
- **[volesti]** `src/volesti/include/volume/volume_cooling_gaussians.hpp` and
  `random_walks/gaussian_ball_walk.hpp`, which is the code ttc runs.
- **[ttc]** `src/volume/lra_volume.cpp` (the oracle call is at line 1255).

## What ttc actually runs
```
volume_cooling_gaussians<GaussianBallWalk, mt19937 seed 3>(P, error = ε/12, walk_length = 10+n²)
```
- `P` is cdd-canonicalised and given the certified Chebyshev ball. No affine transform or
  rounding is applied first.
- The result is clipped to `[(1-ε')(2r/√n)ⁿ, (1+ε')(2R)ⁿ]`.

volesti's code is a C++ port of the [MATLAB] code. The table below shows that it follows the
[MATLAB] almost line for line. It does have two families of absolute constants (rows 2 and 3)
that the MATLAB does not have. The real gap is ttc's use of it: ttc skips the rounding step,
the step that [CV-PVA] and [CV-GC] both require.

## Verdict
For **well-rounded bodies of moderate size** (an inner ball of radius ~0.1–10⁴ and R/r ≲ √n), ttc is
the [CV-PVA] algorithm with three changes: a ball walk instead of hit-and-run, 1 chain instead
of 5, and a larger window. It is accurate there: cubes in 2–10D are within 3–10%, and simplices
up to 8D within 13%.

Outside that regime it is **not** the published algorithm and fails badly:
- **Elongated bodies:** box [0,2]×[0,1000] gives 7.78 instead of 2000. A 5D box with aspect
  ratio 10 gives 0.36× the true volume.
- **Bodies with side ≤ 10⁻³:** the estimate is too low. At side ≤ 10⁻⁴ the oracle never returns
  (infinite loop).
- **Bodies with side ≳ 10⁵:** the estimate is too low. At side 10⁶ the oracle returns `inf` and
  ttc aborts.

Even where it matches the theory, the oracle is a heuristic. The `(1±ε')` with probability
`1-δ'` guarantee that the union-algorithm proof assumes is not delivered (row 7).

## Deviations, ranked by impact

| # | Theory / reference | volesti / ttc | Effect | Status |
|---|---|---|---|---|
| 1 | **Round first.** Step 1 of [CV-PVA] Fig. 2 is `T = Round(K, 8n³)`, then use `vol = \|T\|·…` (§3.1, Fig. 3: SVD of hit-and-run points, repeated until 2-isotropic). The §3.1 motivating example is a 10⁹×1×…×1 box. [CV-GC] assumes `B_n ⊆ K ⊆ C√n B_n` (tex 108, 210). In the [MATLAB] code, `preprocess.m:37` gives the usage example `P = makeBody('long_box',10); P = preprocess(P); vol = Volume(P);` (MVE rounding) | **No rounding anywhere.** volesti has rounding routines (`preprocess/*rounding*`), but ttc never calls them | Mixing needs ~n²R²/r² steps ([CV-PVA] §2, Thm 1). On a skinny body the walk never leaves the centre. The schedule then "sees" a flat Gaussian and stops early (row 4), so the last ratio is estimated on a tiny neighbourhood. Result: a large **under**estimate that grows with R/r | **measured**: 35/67 e2e cases fail, and the long-box failures are all of this kind |
| 2 | **Starting Gaussian search uses a relative tolerance.** [MATLAB] `ConvexBody.m:186` has `while upper-lower > upper*1e-7` | `volume_cooling_gaussians.hpp:67,100` has an **absolute** `tol = 1e-7` | a₀ ∝ 1/r². When r ≲ 10⁻⁴, a₀ ≳ 10⁹, and the gap between adjacent doubles there exceeds 1e-7, so the bisection **never terminates** | **measured**: side 1e-6 hangs (>300 s) |
| 3 | **Ball-walk step ∝ min(σ, r)/√n.** [CV-GC] Fig. "Volume": δ = min{σ,1}/(…√n) with the unit ball inside K, i.e. in units of r | `gaussian_ball_walk.hpp:52` has δ = 4r/√(max(1,a)·n). The `1` is absolute, so the formula is scale-correct only when r = 1. The scale-free form is 4r/√(max(1, a r²)·n) = 4·min(r, 1/√a)/√n | Tiny bodies: δ ∝ r² instead of r, so steps are far too short. Huge bodies: δ ≫ σ, nearly every proposal is rejected, the first schedule check fails, and `a_vals[0]` is replaced by 0, giving `(π/0)^{n/2} = inf` (line 261 together with line 362) | side 1e6 → `inf` **measured**; mechanism **hypothesised**, being verified |
| 3b | **Safe annealing step.** [CV-PVA] Lemma 3 shows aᵢ₊₁ = aᵢ(1−1/n) keeps Var/E² < 1 **only for n ≥ 4**. [CV-GC] uses fixed rates 1+1/√n and 1+σ/(C√n) | `get_next_gaussian` (`:153-176`) keeps doubling k while the *sample* relvar < C. The samples are near the centre, so the test is blind to far tails. In 2D, 1−1/n = 0.5 is not safe at all | Long box: the last phase (a 0.517→0) is estimated at 0.0041× its exact value, and the jumps have exact relvar 378–∞. This is the place where 2000 → 7.78 is lost. Capping the step at a′ ≥ a(1−α_max(n,C)), from the log-concavity of aⁿZ(a), gives 7.78 → 842 | **measured** (annealing agent) |
| 4 | **Schedule end.** [CV-GC] cools until σ² > C²n, a bound set by R, then jumps to uniform | [MATLAB]/volesti stop when a 150/ε-step walk sees E[f_next/f_cur] ≤ 1.001 (`:251`), then *replace* the current aᵢ by 0 (`:261`). This matches `ConvexBody.m:806-811` | The test depends on the data. If the chain hasn't explored K, it declares flatness too early. That is harmless after rounding and is the trigger for row 1's failure without it | matches [MATLAB]; fragile without rounding |
| 5 | **Walk.** [CV-PVA] uses hit-and-run with exact chord sampling (§3.4), and the [MATLAB] default is coordinate hit-and-run (`walk_type=0`). Every step is a sample, with 5 round-robin chains (§3.5) | ttc uses the Gaussian **ball walk** with fixed δ. Each sample takes 10+n² walk steps, and there is 1 chain | Hit-and-run crosses a whole chord in one step. The ball walk moves at most ~4r/√n per step, so it is far more sensitive to R/r. With 1 chain the samples are more dependent | design choice |
| 6 | Sliding window W = 4n²+500 (Fig. 2, `Volume.m:249`) | W = 6n²+800 (`:275`) | More conservative, so convergence is slower but more trustworthy | benign |
| 7 | **Accuracy is heuristic.** [CV-PVA] §3.3: "we can have false positives". In its §4.2 tables every mean is **below** the true volume (2–20% low), and std/mean exceeds ε on simplices (0.31–0.69 against ε = 0.2). [CV-PVA] §4.2.3 recommends averaging several runs | ttc treats the output as a (1±ε/12) oracle with probability 1-δ' (union proof). The oracle seed is fixed (3), so repeat runs cannot be averaged, and δ' is never used | The ε' target is not guaranteed. Expect a systematic low bias of a few %. The clipping bounds hide failures: in the e2e suite, 2–10× underestimates were lifted to the lower clip | known; documented |
| 8 | Per-phase error ε/√m with ε := (1−frac)ε after choosing a₀ (`ConvexBody.m:205`) | `curr_eps = error/√mm` uses the full `error` (`:383`) | 10% looser per-phase target | minor |

## What matches the reference exactly
- a₀ from the Lemma 2 union bound with `frac = 0.1` (`:59-117`, compare [CV-PVA] §3.2.1). This holds
  apart from the tolerance in row 2. The [MATLAB] code also has an ellipsoid term that doesn't
  apply to a polytope.
- Annealing step `a_{i+1} = aᵢ(1−1/n)^k` with doubling and halving of k, the variance test
  `Var/E² ≥ C` with **C = 2** (`:165`). The paper's eq. (1) says < 1, but [MATLAB] uses 2. The
  test uses N = 500C + n²/2 samples, plus a `1+1e-5` progress test.
- The schedule check uses ⌈150/((1−frac)ε)⌉ steps with tolerance 1.001 (`:214,251`).
- The ratio estimator is a running mean of f_{i+1}/f_i along the chain (Fig. 2). The volume is
  `(π/a₀)^{n/2} ∏ Rᵢ` (`:362`). There is no |T| factor because there is no rounding.
- The sliding-window rule `(max−min)/max ≤ ε'/2` (`:430`, `Volume.m:296`).

## Implications for the fix (to be confirmed by the test agents)
1. **Normalise scale** so that the inner radius is 1: x ↦ x/r, then multiply the volume by
   rⁿ, accumulated in log space. This makes rows 2 and 3 harmless, and is what [CV-GC] assumes
   (`B_n ⊆ K`). Patching volesti so it is scale-free (`upper*1e-7` and δ = 4·min(r, 1/√a)/√n) is
   the belt-and-braces version.
2. **Round** before the oracle, then multiply by |det T|. Use either volesti's SVD/isotropic
   rounding, which matches [CV-PVA] Fig. 3, or MVE rounding, which matches [MATLAB]
   `preprocess.m`. Iterate until the singular-value ratio is ≤ ~2, and keep the certified
   inner ball of the *transformed* body.
3. Keep the clipping as a safety net, but log every time it is hit. Consider averaging or
   re-seeding per polytope, following [CV-PVA] §4.2.3.

Per-agent evidence, with exact numbers, is collected in `agent-todo.md`.
