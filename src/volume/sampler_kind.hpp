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
// of the rounded body): 10 + dim^2, but at least 30 for the ball walk, whose
// step 4r/sqrt(n) is too short in low dimension (2-D cube needs ~28 steps;
// tests/volesti/test_sampling_mixing.cpp).
inline unsigned defaultSampleWalkLength(SamplerWalk walk, std::size_t dimension)
{
  const unsigned base = static_cast<unsigned>(10 + dimension * dimension);
  return walk == SamplerWalk::Ball ? std::max(base, 30u) : base;
}

}  // namespace ttc
