#pragma once

// The sampling walk selectable with --sampler.  Kept free of volesti includes
// so lra_volume.hpp / main.cpp stay light; the sampler itself lives in
// polytope_sampler.hpp.

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

}  // namespace ttc
