# Union accuracy test plan, by component

*2026-11-14. Source: `text-for-agents/kr25-ssm` (main-smt-vol.tex, sections/algo.tex,
analysis.tex, proof-appendix.tex). The single-polytope pieces (volume oracle `voled`,
`GenerateSamples`) are already tested (`tests/volesti/`, see `agent-todo.md`). This plan covers
everything else that decides whether the **union** estimate is right.*

## Theorem → components

Theorem `volumettc`: with probability at least 1−δ, est ∈ (1±ε)·vol(F). The proof is a chain.
Each link is one component, and each component gets its own tests:

| # | Component (Alg. 1 line) | Code | Lemma it must satisfy | What breaks if it is wrong |
|---|---|---|---|---|
| C1 | Boolean abstraction + AllSAT → cubes → polytopes (l.1–2) | `allsat_volume.cpp`, `buildPolytopeMatrices`, `canonicalizePolytope` (cdd) | `thm:decomposition`: ⋃ P(cube_i) = sol(F); every P_i ⊆ sol(F) | Missing or extra region, so a systematic error |
| C2 | GetPrecision (l.3, Alg. 2) | `precisionFromCubes` | `lem:lattice`, `lem:nonconv`: (1−η)\|Q∩L\|/10^{bn} ≤ vol(Q) ≤ (1+η)\|Q∩L\|/10^{bn}, with η = ε/96 | Lattice count ≠ volume |
| C3 | ε′, δ′, Thresh (l.4–5) | `lraSampleThreshold` | Thresh formula; ε′ = ε/12; δ′ = δ/(4m) | Too few samples, so failure rate above δ |
| C4 | Volume oracle + clipping (l.7, clip to [(1−ε′)(2γ/√n)ⁿ, (1+ε′)(2R)ⁿ]) | `RoundedPolytope`, `certifiedInnerBall`, `outerRadius` | voled ∈ (1±ε′)vol, and the clip bounds are sound | Wrong λ in the Poisson invariant |
| C5 | Remove covered points (l.8–10) | `SampleStore::removeInside` (double / GMP) | Exact membership in P_i | Double counting of overlaps |
| C6 | Halving (l.11–14) and the threshold loop (l.16–19), Poisson/Binomial | `runUnionAlgorithm`, `downsampleHalf` | `lem:process` (three equivalent procedures); invariant k_x ~ Pois(λ_j) | Biased \|X\|, wrong variance |
| C7 | GenerateSamples (l.20) | `sampleLattice`, `sampleRestartWalk` | eq. (almostuniform), plus independence | Overlap regions mis-weighted |
| C8 | Output (l.22) and the error events E_j, Noise_j | `result.volumeEstimate` | `lem:numlattice`: est ∈ (1±ε/2)\|Q∩L\|/10^{bn} with probability ≥ 1−δ | Everything above, compounded |

C4 and C7 are covered by the volume and sampling work (`ttc-volume-vs-theory.md`,
`ttc-sampling-vs-theory.md`). The rest are new.

## Prerequisite refactor (bit-identical)
`runUnionAlgorithm`, `SampleStore`, `PolytopeVolume`, `precisionFromCubes` and
`lraSampleThreshold` move from `lra_volume.cpp` into `src/volume/union_estimator.{hpp,cpp}`.
- The union loop stays a header template on `Coord` and on the sample generator, so tests can plug
  in **exact** components.
- This is the same pattern as `polytope_sampler.{hpp,cpp}`.
- Acceptance: the 16-run `s volume` bit-identity check used for the earlier refactors.

## Tests

### T1. Union estimator with exact components (C6, C8): `test_union_core.cpp`
The heart of the plan. Drive the extracted `runUnionAlgorithm` with:
- **exact volumes** (boxes, simplices, polygons);
- an **exact uniform lattice sampler**: enumerate P∩L for small b, or draw uniform continuous
  points and lattice-round them.

With every other error source removed, the estimator alone must meet Lemma `numlattice`.

| Case | Exact reference | Stresses |
|---|---|---|
| 2 identical boxes | \|B∩L\| | removal of duplicates |
| Nested boxes A ⊂ B, both orders | \|B∩L\| | order dependence, removeInside |
| k disjoint boxes, k = 2…50 | Σ\|B_i∩L\| | many polytopes; Thresh uses m |
| Chain of k overlapping boxes | inclusion–exclusion | overlap bookkeeping |
| Volumes 1e-6 … 1e6 in one union | Σ | many halvings, log2 p underflow |
| m = 200 tiny + 1 huge, and reversed order | exact | `while p ≥ Thresh/t` path |
| Thresh boundary: N_i + \|X\| ≈ Thresh | exact | Binomial thinning loop |

- **Accuracy:** for each case, ≥ 200 seeds at (ε, δ) = (0.8, 0.2) and (0.3, 0.1). The fraction of
  runs with |est/exact − 1| > ε/2 must be ≤ δ (one-sided binomial test, α = 0.01).
- **Unbiasedness:** |mean(est)/exact − 1| ≤ 3·SE.
- **Distribution:** the store's per-point copy counts follow Pois(λ_j) (dispersion index ≈ 1), as
  in `lem:process`.
- **Order invariance:** the estimate distributions of permuted polytope orders agree (two-sample KS).

### T2. Error injection (C4, C7 budget): `test_union_robust.cpp`
Same driver, with controlled imperfections:
- each volume ×(1±ε′), fixed or random sign;
- a sampler with total-variation distance η: mix the uniform sampler with a point mass, or a
  centre-biased law.

Expected: est stays within (1±ε) while volume errors are ≤ ε′ and sampler TV ≤ ε′/m. The bias grows
linearly beyond that. This measures the real tolerance of the MCMC-based sampler; ttc-sampling
row 3 notes the lemma's pointwise assumption is not met.

### T3. Lattice counting vs volume (C2): `test_union_lattice.cpp`
- Recompute `precisionFromCubes` by hand from Alg. 2: γ = 64mn/η·√log(4r/η), b = ⌈log₁₀(γ/γ̂_min)⌉.
- For small 2D/3D unions (b ≤ 3), enumerate |Q∩L| exactly and check `lem:lattice`:
  |10^{−bn}|Q∩L| / vol(Q) − 1| ≤ η.
- Also probe the lemma's slack: the smallest b that already gives η-accuracy, against b from
  GetPrecision. How conservative is it, and what does it cost (b drives GMP digits and C)?

### T4. Parameters (C3, C4 clipping): `test_union_params.cpp`
- `lraSampleThreshold` equals the l.5 formula for a grid of (ε, δ, m).
- ε′ = ε/12 and δ′ = δ/(4m) as used.
- Clipping is sound: (1−ε′)(2γ/√n)ⁿ ≤ vol ≤ (1+ε′)(2R)ⁿ on box, simplex, cross-polytope,
  rotated and skinny bodies. outerRadius R ≥ the true circumradius about the centre.

### T5. Membership and the store (C5, C6 mechanics): inside `test_union_core.cpp`
- `removeInside` is exact on points on facets and on lattice points of shared faces, in double
  and GMP. Use the integer-units rule from the sampler fix.
- `downsampleHalf` keeps each point independently with probability ½ (binomial test).
- The Poisson draw has mean t·p (large-λ regime).
- `log2Density` bookkeeping gives the same result as exact p·C computed in GMP.

### T6. Decomposition (C1): `tests/volesti/union_smt2/` + `run_union_decomp.py`
- With `--dump-ine PREFIX`, read the dumped polytopes. For random points in a bounding box (and in
  slabs near each facet), check that F(x) ⇔ ∃i: x ∈ P_i, evaluating F with a small SMT-LIB
  evaluator in Python.
- Each P_i must be ⊆ sol(F) (sample inside P_i, evaluate F).
- Formulas exercise `or`/`and`/`not`/`ite`, `=` (zero-volume pieces), strict and non-strict
  inequalities, and redundant literals.
- Report cover size and overlap ratio Σvol(P_i)/vol(Q). Overlap is what C5/C6 must handle.

### T7. End-to-end union accuracy (all): `tests/volesti/union_smt2/*.smt2` + `run_smt2.py --seeds K`
QF_LRA formulas with **exactly known** vol(sol F), each with an `; exact-volume:` header:

| Family | Exact volume by |
|---|---|
| Union of 2–10 overlapping boxes | inclusion–exclusion |
| Nested / identical / disjoint boxes | direct |
| L-shape, frame (square minus square), staircase | direct |
| Union of two triangles forming a square; hexagon split into 6 triangles | direct |
| Skinny overlapping slabs (aspect 100–1000), rotated | inclusion–exclusion / det |
| Many tiny disjoint boxes + one large box | Σ |
| Dims 2–8; sides 1e-3 … 1e3 | scaling |
| Parity-like Boolean structure (2^{v−1} cubes, `thm:decompose`) | combinatorial |

- **Per formula:** K = 20 seeds. Pass if ≤ ⌈δK⌉ + binomial slack runs fall outside (1±ε), and
  the mean is within ε/2.
- **Run with:** default; `--nogmp`; `--fullgmp`; `--sampler cdhr`.
- **Regression:** run the same suite with the start binary (`update` @ `75e05e0`) to show what the
  fixes changed.

## Order of work and cost
1. Refactor (bit-identical): ~1 h.
2. T1 + T5: the core estimator. These decide whether the union logic itself is correct.
3. T7: end-to-end, old vs new.
4. T3, T4: cheap, deterministic.
5. T2: tolerance study, informs the paper's lemma.
6. T6: needs a small SMT-LIB evaluator.

T1, T2 and T7 run many seeds and are CPU-heavy, so give them `slow` tags.
