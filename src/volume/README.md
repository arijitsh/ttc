# LRA volume engine

ttc estimates the volume of the solutions of an LRA formula F with Algorithm 1
of the manuscript (`kr25-ssm/sections/algo.tex`): split F into polytopes, then
estimate the volume of their union with a streaming sampler. The driver,
`unionvolume.cpp`, is short and names each step. This table maps every line
of the pseudocode to the file that implements it and the unit test that
covers it.

| Alg. 1 | Pseudocode | Code | Tests (`tests/lra/`) |
|---|---|---|---|
| 1 | `F_bc, M <- BoolAbs(F)` | `boolean_abstraction.hpp` (`buildBooleanAigFromTerm`), called by `TTCParser::computeCubes` | `test_boolabs.cpp` |
| 2 | `F_bd <- ToDNF(F_bc)` | `dnf_enumeration.hpp` (`enumerateCubes`), using AllSAT from `aig_allsat.hpp` | `test_todnf.cpp` |
| | `P(c_i)` for each cube | `linear_constraints.hpp` (literals to `A x <= b`), `canonicalize.hpp` (cdd), `prepared_polytope.hpp` (balls, drop zero volume) | `test_cube_polytope.cpp`, `test_polytope_balls.cpp` |
| 3 | `b <- GetPrecision(F_bd, M, eps/96)` (Alg. 2) | `get_precision.hpp` (`getPrecision`, `latticePrecision`) | `test_get_precision.cpp` |
| 4-5 | `eps', delta', Thresh` | `union_parameters.hpp` (`unionParameters`, `lraSampleThreshold`) | `test_union_parameters.cpp` |
| 6-7 | `p <- 1; X <- {}`, loop over `i` | `union_estimator.hpp` (`estimateUnionVolume`, `ScaledProbability`) | `test_union_estimator.cpp` |
| 8 | `t <- Thresh 10^{bn} VolEst(P_i, eps', delta')` | `vol_est.hpp` (`VolEst`: oracle, retries, clipping), with `polytope_volume.hpp` (rounding and VolEsti oracles) | `test_vol_est.cpp` |
| 9-11 | remove `s in X` with `s in P_i` | `sample_multiset.hpp` (`removeInside`) | `test_sample_multiset.cpp` |
| 12-18 | halve `X` and `p`; `N_i <- Poisson(t p)`; thin `N_i` | `union_estimator.hpp`, `sample_multiset.hpp` (`halve`) | `test_sample_multiset.cpp`, `test_union_estimator.cpp` |
| 19 | `S <- GenerateSamples(P_i, N_i, b, eps', delta')` | `generate_samples.hpp` (inflate, propose, round, reject), with proposals from `polytope_sampler.hpp` | `test_generate_samples.cpp` |
| 20, 22 | `X.Append(S)`; output `\|X\| / (p Thresh 10^{bn})` | `sample_multiset.hpp`, `union_estimator.hpp` | `test_sample_multiset.cpp` |
| all | the whole algorithm | `unionvolume.hpp` (`computeLraVolume`) | `test_unionvolume.cpp` |

Supporting files:

- `cube.hpp`: a cube `c_i`, as the theory literals it fixes.
- `polytope_balls.hpp`: the certified inner ball and the outer radius,
  `B(z, gamma) ⊆ P ⊆ B(z, R)`, computed with LPs in lp_solve.
- `volume_options.hpp`: `VolumeOptions` (the CLI knobs), `VolumeSetup`, and the
  choice between double and GMP arithmetic (`GmpMode`).
- `seed_streams.hpp`: per-polytope random streams derived from `--seed`.
- `mp_float.hpp` and `gmp_walks.hpp`: GMP floats and the walks used with `--fullgmp`.
- `ine_writer.hpp`: `--dump-ine`.
- `volesti_body.hpp`: internal to the `.cpp` files. It is the only place
  volesti's polytope type appears; no public header includes volesti.

## Build and test

`CMakeLists.txt` builds these sources once, as the OBJECT library
`ttc_lra_volume`. `ttc` and the doctest binary `lra_unit_tests` both link it:

```
make -C build-gcc -j16 ttc lra_unit_tests
build-gcc/lra_unit_tests            # about 10 s; ctest runs it as `lra_unit`
```

`tests/volesti/` holds the longer statistical audits of the same components
(`make -C tests/volesti`).

## Where the code departs from the paper

- **VolEst** runs a practical oracle (VolEsti's Cousins-Vempala cooling, by
  default with coordinate hit-and-run). It has no confidence parameter, so
  `delta'` is unused. Failed runs are retried on fresh streams. The result is
  clipped to `[(1-eps')(2 gamma_i/sqrt n)^n, (1+eps')(2 R_i)^n]`.
- **GenerateSamples** proposes points with a billiard walk of `n` steps,
  restarted at the Chebyshev centre of the rounded body (`--sampler` and
  `--walklen-samp` change this). This is a practical walk length, not a
  certified mixing time. It samples in double while
  `b + log10 max|x| <= 15`, otherwise in GMP (`--gmp`, `--fullgmp`, `--nogmp`).
- **GetPrecision** takes `r` and `min gamma_hat` over the polytopes of positive
  volume only; `m` counts every cube.
- **Line 8** runs inside the loop, as in the pseudocode, so each table row is
  printed as soon as its polytope is done. VolEst draws from its own streams,
  so where it runs does not change any estimate.
- **BoolAbs** uses cvc5's in-memory AIG (`Solver::getBooleanAbstractionAig`)
  when the linked cvc5 provides it. The cvc5 pinned in `configure.sh` does
  not, so ttc builds the AIG itself (`buildBooleanAigFromTerm`).
- **Known bug:** if even the inscribed ball of a cube's polytope is unbounded
  (for example `x >= 0 & y >= 0`), the Chebyshev LP returns UNBOUNDED, the cube
  is dropped as flat, and ttc prints `s volume 0` instead of an error. The
  `should_fail` case in `test_cube_polytope.cpp` tracks this.
