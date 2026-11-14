# How close is ttc's single-polytope sampler to what the theory needs?

*2026-11-14. Scope: drawing samples from **one** polytope P = {x : Ax ≤ b} (GenerateSamples), not the
union estimator. Companion to `ttc-volume-vs-theory.md`. Evidence: the s-theory agent's audit
(notes and simulation in `scratchpad/s-theory/`), `tests/volesti/test_sampling_smoke.cpp`, and the
sampling agents' tests as they report (`agent-todo.md`, step 12).*

## Sources
- **[SM]** ttc manuscript, camera-ready copy `~/Dropbox/research/tex/SAT25-SM_cam/sections`.
  - GenerateSamples is in `algo.tex`.
  - The lattice lemma (lem:nonconv) is in `analysis.tex`.
  - The Poisson-splitting proof is in the proof appendix, l.17–39.
  - `implementation.tex` l.30–34 covers hit-and-run.
  - The revised constants the code cites (η = ε/96, "cor:linf" κ) are *not* in this copy; they were
    checked against the code comments.
- **[CV]** Cousins & Vempala's MATLAB sampler (`Volume-and-Sampling-master.zip`):
  - `polySampler.m` and `Sample.m` round first (`preprocess.m`), burn in, then keep **one chain**
    with thinning (8n² coordinate hit-and-run steps between samples).
- **[PG]** Polyak & Gryazina 2014, the billiard walk. It proves stationarity of the uniform
  distribution only; there is **no mixing bound**.
- **[LV]** Lovász–Vempala hit-and-run.
  - From a warm start it needs O*(n²R²/r²) steps.
  - From a fixed point it needs O*(n³R²/r²), with constants above 10¹⁰.
- **[volesti]** `random_walks/uniform_*_walk.hpp`, `compute_diameter.hpp` and
  `sampling/random_point_generators.hpp`.

## What ttc does (`src/volume/polytope_sampler.{hpp,cpp}`)
1. `sampleRestartWalk`: for **every** sample, restart the chosen walk at the certified Chebyshev
   centre c, take walkLength = 10+n² steps, and keep the endpoint. The polytope is **not rounded**.
   The default walk is billiard, with trajectory length L = 2√n·r (`compute_diameter.hpp:39`).
2. `sampleLattice`, the GenerateSamples step:
   - Inflate each facet by 10⁻ᵇ·κ·‖A_j‖, where κ = √(2 ln(16/ε′)) + √(2 ln m).
   - Draw proposals from the inflated body.
   - Round each coordinate stochastically to 10⁻ᵇℤ.
   - Reject points outside P, with a retry budget derived from δ.

## What the proof needs
- **Pointwise near-uniformity on the lattice:**
  Pr(s) ∈ [(1−η)/|S|, (1+η)/|S|] for every s ∈ P ∩ 10⁻ᵇℤⁿ.
- **Independence across draws**, which the Poisson-splitting argument relies on.
- The lattice lemma also assumes the **continuous proposal is exactly uniform** on the inflated
  body.

## Verdict
- **Independence: met.** Restarts make the draws i.i.d.
- **Distribution: not met** once the body is skinny. Each draw follows the
  (walkLength+1)-step distribution from c, which can never leave a ball of radius about
  (10+n²)·2√n·r around c.
  - Long box [0,2]×[0,1000]: **100% of the samples land in y ∈ [900,1000]**, against 10% for a
    uniform sampler. All of them lie within about ±17 of the centre y=999.
  - Simulation: KS distance from uniform is 0.41 (n=5) and 0.36 (n=10) at aspect ratio 100;
    aspect 10 is at noise level.
- **Lattice step: matches** the manuscript's construction and constants.
- **The pointwise-uniformity assumption cannot be met** by any MCMC sampler. Every MCMC
  guarantee, even an asymptotic one, is about total variation.

## Deviations, ranked by impact

| # | Theory / reference | ttc | Effect | Status |
|---|---|---|---|---|
| 1 | Samples (near-)uniform on P. [CV] and volesti's own tools keep one rounded, burnt-in chain with thinning | Restart at a fixed centre with a fixed 10+n² steps (`polytope_sampler.cpp`, `sampleRestartWalk`) | Each sample is the short-walk law from c, not uniform. The bias is heavy once R/r ≳ 30. In the union algorithm, overlaps between polytopes are mis-estimated (up to 100% when one polytope's samples never reach the overlap) | **measured** (smoke test, simulation) |
| 2 | Round first ([CV] `polySampler.m`; every mixing bound scales with (R/r)²) | Not rounded, although the volume oracle now is (`polytope_volume.cpp`, `RoundedPolytope`) | Same root cause as the volume bug: steps scale with the inner radius r | **measured** |
| 3 | [SM] lem:nonconv: pointwise (1±η) on lattice points, with the proposal exactly uniform | MCMC gives at best total-variation closeness | The proof does not cover the implementation even with perfect mixing. **Paper-side fix:** restate the lemma in TV terms. A TV error ε_TV per polytope shifts the expected count by at most ε_TV·N_i, so ε_TV ≤ ε′/m suffices for the expectation; a coupling gives approximate independence | theory gap |
| 4 | Walk length from a mixing bound or a diagnostic ([CV]: burn-in, then 8n² steps between samples; ESS/PSRF or a sliding window) | Fixed 10+n², with no diagnostic. `sampleLattice`'s retry budget cannot detect bias, because a stuck walk almost never proposes outside P | No safety net | design |
| 5 | Lattice cells resolved exactly | `--nogmp` and the default point-repr mode run the walk in double. Once ulp(\|x\|) approaches 10⁻ᵇ (e.g. \|x\| ≈ 10³ with b ≥ 12), proposals sit on a grid as coarse as a lattice cell | Cell probabilities stop being uniform. The only guard is `spacing == 0` | **hypothesised** |
| 6 | GMP walks (`--fullgmp`) | Same restart kernels. Rows of A are not normalised, and the facet test uses an absolute `1e-12`, so rows with tiny norms can be ignored. Accelerated billiard falls back to billiard. Each sample takes t steps, not volesti's t+1 | GMP fixes arithmetic, not mixing. Issue 1 applies unchanged | **source-checked** |
| 7 | Billiard reversibility | volesti and the GMP walks both give up after 50n reflections and keep the start point (`uniform_billiard_walk.hpp:107,126`) | Slight non-reversibility; rare in round bodies | minor |

## What matches
- **Lattice step.** It uses ℓ₂ row norms and the same 10⁻ᵇ for inflation and the rounding cell.
  Each coordinate gets Kannan–Vempala randomized rounding (Bernoulli of the fractional part), and
  the rounding RNG is seeded separately from the walk. Rejection is against the original b.
- **Per-polytope facet count.** Using r_i facets in κ is valid for the per-polytope union bound.
- **κ is large enough.** It uses ln(16/ε′) = ln(192/ε), slightly below the ln(4/η) = ln(384/ε)
  that η = ε/96 suggests (3.31 vs 3.51 at ε = 0.8). It still suffices because Hoeffding gives
  exp(−2κ²).
- **Retry budget.** The bound n·(ε′/3)ᵏ ≤ δ is a correct union bound, assuming uniform proposals.

## What was changed (2026-11-14)
| Deviation | Fix | Commit |
|---|---|---|
| 1, 2. Restart on the unrounded body | The walk runs in the rounded body, using `roundPolytope`, the same max-inscribed-ellipsoid map as the volume oracle, including in `--fullgmp`. Each point is mapped back with x = shift + T·y before the lattice step. The restart design stays: the mixing tests showed that one long chain keeps lag-1 correlation | `08026a2` |
| 4. Walk length | Per-walk defaults (`sampler_kind.hpp`): billiard 10+n²; ball walk ≥ 100; RDHR/CDHR 8·(10+n²), all chosen from the uniformity tests | `08026a2`, `722ff7b` |
| new | Billiard bias: volesti reflects at `dl = 0.995` of the distance to the facet, which under-samples a boundary layer at every walk length. Changed to 0.999999 | `722ff7b` |
| new | Lattice membership: `A(k/10ᵇ) ≤ b` dropped lattice points lying exactly on a facet (floating-point noise). Now tested as `A·k ≤ 10ᵇ·b`, exact in GMP | `08026a2` |
| new | A tiny body far from the origin gave samples *outside* P (slack cancellation). Fixed by rounding, which centres the body | `08026a2` |

Still open:
- Row 3: paper-side; the lattice lemma should be restated with a TV hypothesis.
- Row 5: no double-precision guard yet.
- Row 6: GMP walk tolerances.
- Row 7: the 50n-reflection fallback.

## Recommended fix (written before the fix)
1. **Round before sampling.** Reuse the max-inscribed-ellipsoid rounding that the volume oracle
   already uses (`RoundedPolytope`): walk in y-space, then map back with x = T·y + shift *before*
   the lattice step. Inflation, stochastic rounding and rejection stay in x-space. Affine maps
   preserve uniformity.
2. **One chain per polytope.** Burn in from the centre, then thin by k steps between samples.
   After rounding, k = 10+n² billiard steps is the natural starting point; check it with the
   mixing agent's gauge-KS numbers.
3. **Precision guard.** Require max|x|·10ᵇ·2⁻⁵² ≪ ε′/n, or switch to the full GMP walk.
4. **Paper.** Restate the lattice lemma with a total-variation hypothesis (row 3).
