#pragma once

// The sampling walk selectable with --sampler.  Kept free of volesti includes
// so lra_volume.hpp / main.cpp stay light; the sampler itself lives in
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
//   billiard / accelerated billiard: 10 + dim^2 (passes from 1-2 steps once
//     the body is rounded; test_sampling_mixing.cpp);
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
    default: return base;
  }
}

}  // namespace ttc
