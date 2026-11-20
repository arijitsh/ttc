#pragma once

// The sampling walk selectable with --sampler.  Kept free of volesti includes
// so volume_options.hpp / main.cpp stay light; the sampler itself lives in
// polytope_sampler.hpp.

#include <algorithm>
#include <cstddef>

namespace ttc
{

// Random-walk used by volesti to draw (approximately) uniform points from each
// polytope during the union-volume estimation. Selected with --sampler.
//   Billiard (default):    uniform_billiard_walk.hpp -- reflecting trajectories,
//                          mixes fastest for H-polytopes.
//   AcceleratedBilliard:   uniform_accelerated_billiard_walk.hpp.
//   Ball:                  uniform_ball_walk.hpp -- classic ball walk.
//   RDHR / CDHR:           random- / coordinate-directions hit-and-run.
// Under --fullgmp, accelerated billiard falls back to GMP billiard.
enum class SamplerWalk
{
  Billiard,
  AcceleratedBilliard,
  Ball,
  RDHR,
  CDHR
};


// Default number of walk steps per sample (each sample restarts at the centre
// of the rounded body), from the uniformity tests in tests/volesti/:
//   billiard / accelerated billiard: dim (passes from 1-2 steps once the body
//     is rounded; test_sampling_mixing.cpp).  No mixing bound is known for
//     the billiard walk; dim steps is at least the O(sqrt(dim)) warm-start
//     burn-in of Chalkis et al. (SoCG 2021) and the single step VolEsti uses
//     (Chalkis-Fisikopoulos 2021).  On krbenchsmt, 10 + dim^2 steps changed
//     no estimate beyond run-to-run noise while sampling 8-9x slower;
//   ball walk: at least 100 -- its step 4r/sqrt(n) overshoots a rounded
//     low-dimensional body, most proposals are rejected and a short restart
//     can stay at the centre (atom; test_sampling_lattice.cpp failed at 30);
//   hit-and-run (RDHR, CDHR): 8 (10 + dim^2) -- rounded simplices' and
//     cross-polytopes' sharp vertices are under/over-covered at 4x
//     (test_sampling_shapes.cpp).
inline unsigned defaultSampleWalkLength(SamplerWalk walk, std::size_t dimension)
{
  const unsigned base = static_cast<unsigned>(10 + dimension * dimension);
  switch (walk)
  {
    case SamplerWalk::Ball: return std::max(base, 100u);
    case SamplerWalk::RDHR:
    case SamplerWalk::CDHR: return 8 * base;
    case SamplerWalk::Billiard:
    case SamplerWalk::AcceleratedBilliard:
    default: return static_cast<unsigned>(std::max<std::size_t>(1, dimension));
  }
}

}  // namespace ttc
