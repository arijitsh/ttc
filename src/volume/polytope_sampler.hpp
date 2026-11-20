#pragma once

// The continuous sampler inside GenerateSamples (generate_samples.hpp): points
// of a polytope from a random walk restarted at its Chebyshev centre for each
// point.  Restarting, rather than thinning one chain, keeps the points
// independent, as the analysis requires.
//
// With `round` (ttc's default) the walk runs in the rounded body y
// (roundPolytope: maximum inscribed ellipsoid -> unit ball) and each point is
// mapped back, x = shift + T y; affine maps preserve uniformity.  Unrounded,
// every walk steps in units of the inner radius, so on a skinny body
// (R/r >> sqrt n) the samples never leave the centre (long box
// [0,2]x[0,1000]: 100% of samples in the last tenth), and a tiny body far
// from the origin loses its slacks to cancellation (samples outside P).

#include <vector>

#include <Eigen/Dense>

#include "volume/mp_float.hpp"
#include "volume/sampler_kind.hpp"

namespace ttc
{

// n points of {x : A x <= b}, each the endpoint of a fresh `walkLength`-step
// run of volesti's `walk` from the certified Chebyshev centre, on an mt19937
// stream seeded with `seed`.  ballSeconds receives the time spent on the
// inner ball and the rounding.
std::vector<Eigen::VectorXd> sampleRestartWalk(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                               long long n, SamplerWalk walk, unsigned walkLength,
                                               unsigned seed, double& ballSeconds,
                                               bool round = true);

// The same in GMP arithmetic (--fullgmp), with the walks of gmp_walks.hpp;
// accelerated billiard has no GMP version and runs the GMP billiard walk.
std::vector<MpVector> sampleRestartWalkGmp(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                           long long n, SamplerWalk walk, unsigned walkLength,
                                           unsigned seed, double& ballSeconds, bool round = true);

}  // namespace ttc
